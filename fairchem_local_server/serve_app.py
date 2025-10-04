import time

from fastapi import FastAPI
from ray import serve

from .model_runtime import (
    MODEL_NAME,
    TASK_NAME,
    UMA_DEPLOYMENT_NAME,
    _PredictDeploy,
    install_predict_handle,
)
from .models import MDIn, RelaxIn, SimpleIn
from .services import md_step, relax, simple_calculate

app = FastAPI(title="UMA Serve API", debug=True)


@serve.deployment(ray_actor_options={"num_gpus": 0})
@serve.ingress(app)
class Ingress:
    def __init__(self, predict_handle):
        # Inject the UMA handle; builds and caches FAIRChemCalculator.
        install_predict_handle(predict_handle)

    @app.get("/serve/health")
    def health(self):
        return {"status": "ok"}

    # SYNC endpoints so Starlette runs them in a threadpool,
    # which makes the blocking FAIRChem/ASE code safe.
    @app.post("/serve/simple")
    def simple_ep(self, inp: SimpleIn):
        return simple_calculate(inp)

    @app.post("/serve/relax")
    def relax_ep(self, inp: RelaxIn):
        return relax(inp).dict()

    @app.post("/serve/md")
    def md_ep(self, inp: MDIn):
        return md_step(inp).dict()


def deploy():
    # Ray + Serve will be started implicitly by serve.run if needed.
    uma = _PredictDeploy.options(name=UMA_DEPLOYMENT_NAME).bind(MODEL_NAME, TASK_NAME)
    dag = Ingress.bind(uma)
    serve.run(dag, name="http_app", route_prefix="/")
    return dag


if __name__ == "__main__":
    deploy()
    while True:
        time.sleep(60)

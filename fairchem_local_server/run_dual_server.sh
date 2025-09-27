#!/usr/bin/env bash
# run_dual_server.sh - Launch UMA (FAIRChem) HTTP (8000) + HTTPS (8444) servers.
# Generates a self-signed certificate with SANs if none exists.
#
# Usage:
#   ./run_dual_server.sh            # default model/task
#   UMA_MODEL=uma-s-1p1 UMA_TASK=omol ./run_dual_server.sh
#
# Optional env vars:
#   UMA_HTTP_PORT (default 8000)
#   UMA_HTTPS_PORT (default 8444)
#   UMA_CERT_DIR (default ./certs)
#   UMA_CERT_DAYS (default 365)
#   UMA_HOSTNAMES (comma list for SAN, default: kalman,localhost)
#   UMA_IPS (comma list for SAN IPs, default: 127.0.0.1,192.168.1.141)
#   UMA_ONLY_HTTPS=1 (skip HTTP listener)
#   UMA_ONLY_HTTP=1  (skip HTTPS listener)
#   UMA_REGEN_CERT=1 (force regenerate cert)
#
set -euo pipefail

HTTP_PORT=${UMA_HTTP_PORT:-8000}
HTTPS_PORT=${UMA_HTTPS_PORT:-8444}
CERT_DIR=${UMA_CERT_DIR:-$(dirname "$0")/certs}
CERT_DAYS=${UMA_CERT_DAYS:-365}
MODEL=${UMA_MODEL:-uma-s-1p1}
TASK=${UMA_TASK:-omol}
HOSTNAMES_CSV=${UMA_HOSTNAMES:-kalman,localhost}
IPS_CSV=${UMA_IPS:-127.0.0.1,192.168.1.141}
REGEN=${UMA_REGEN_CERT:-0}
ONLY_HTTPS=${UMA_ONLY_HTTPS:-0}
ONLY_HTTP=${UMA_ONLY_HTTP:-0}
KEY_FILE=$CERT_DIR/fairchem-key.pem
CERT_FILE=$CERT_DIR/fairchem-cert.pem
CONF_FILE=$CERT_DIR/san.cnf

mkdir -p "$CERT_DIR"

make_san_conf() {
  local i=1
  echo "[req]" > "$CONF_FILE"
  echo "default_bits = 2048" >> "$CONF_FILE"
  echo "prompt = no" >> "$CONF_FILE"
  echo "default_md = sha256" >> "$CONF_FILE"
  echo "x509_extensions = v3_req" >> "$CONF_FILE"
  echo "distinguished_name = dn" >> "$CONF_FILE"
  echo "[dn]" >> "$CONF_FILE"
  echo "C=US" >> "$CONF_FILE"
  echo "ST=Dev" >> "$CONF_FILE"
  echo "L=Local" >> "$CONF_FILE"
  echo "O=Local Dev" >> "$CONF_FILE"
  echo "CN=${HOSTNAMES_CSV%%,*}" >> "$CONF_FILE"
  echo "[v3_req]" >> "$CONF_FILE"
  echo "subjectAltName = @alt_names" >> "$CONF_FILE"
  echo "[alt_names]" >> "$CONF_FILE"
  IFS=',' read -r -a HARR <<< "$HOSTNAMES_CSV"
  for h in "${HARR[@]}"; do
    echo "DNS.$i = $h" >> "$CONF_FILE"; i=$((i+1))
  done
  IFS=',' read -r -a IARR <<< "$IPS_CSV"
  for ip in "${IARR[@]}"; do
    echo "IP.$i = $ip" >> "$CONF_FILE"; i=$((i+1))
  done
}

if [[ $REGEN -eq 1 || ! -f $KEY_FILE || ! -f $CERT_FILE ]]; then
  echo "[run_dual_server] Generating self-signed certificate (SANs)..."
  make_san_conf
  openssl req -x509 -nodes -days "$CERT_DAYS" -newkey rsa:2048 \
    -keyout "$KEY_FILE" -out "$CERT_FILE" -config "$CONF_FILE" >/dev/null 2>&1
  echo "[run_dual_server] Generated cert: $CERT_FILE"
else
  echo "[run_dual_server] Using existing certificate in $CERT_DIR"
fi

export UMA_MODEL="$MODEL"
export UMA_TASK="$TASK"

PY_MOD=fairchem_local_server.server:app

launch_http() {
  echo "[run_dual_server] Starting HTTP server on :$HTTP_PORT (MODEL=$MODEL TASK=$TASK)"
  uvicorn "$PY_MOD" --host 0.0.0.0 --port "$HTTP_PORT"
}

launch_https() {
  echo "[run_dual_server] Starting HTTPS server on :$HTTPS_PORT (MODEL=$MODEL TASK=$TASK)"
  uvicorn "$PY_MOD" --host 0.0.0.0 --port "$HTTPS_PORT" \
    --ssl-keyfile "$KEY_FILE" --ssl-certfile "$CERT_FILE"
}

if [[ $ONLY_HTTPS -eq 1 && $ONLY_HTTP -eq 1 ]]; then
  echo "[run_dual_server] BOTH UMA_ONLY_HTTP and UMA_ONLY_HTTPS set—choose one." >&2
  exit 1
fi

# Run desired servers. If both, use background + wait.
if [[ $ONLY_HTTP -eq 1 ]]; then
  launch_http
  exit $?
elif [[ $ONLY_HTTPS -eq 1 ]]; then
  launch_https
  exit $?
else
  # dual mode
  launch_http &
  PID_HTTP=$!
  sleep 1
  launch_https &
  PID_HTTPS=$!
  echo "[run_dual_server] HTTP pid=$PID_HTTP HTTPS pid=$PID_HTTPS"
  trap 'echo "[run_dual_server] Stopping"; kill $PID_HTTP $PID_HTTPS 2>/dev/null || true' INT TERM
  wait $PID_HTTP $PID_HTTPS
fi

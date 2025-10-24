import { installBabylonStub } from './helpers/installBabylonStub.js';
import {
  refreshMoleculeMasters,
  getMoleculeMasters,
  getAnyMaster,
} from '../public/vr/vr-utils.js';

describe('x-vr late bond rotation', () => {
  beforeEach(() => {
    installBabylonStub();
  });

  test('masters collapse to molecule_root so late additions inherit rotation', () => {
    const root = new BABYLON.TransformNode('molecule_root');
    root.rotationQuaternion = BABYLON.Quaternion.Identity();
    const earlyMaster = {
      name: 'bond_master_initial',
      getClassName: () => 'Mesh',
      rotationQuaternion: BABYLON.Quaternion.Identity(),
      parent: root,
    };
    root.addChild(earlyMaster);

    const scene = {
      uid: 101,
      meshes: [earlyMaster],
      transformNodes: [root],
    };

    const initial = refreshMoleculeMasters(scene, { force: true });
    expect(initial).toEqual([root]);
    expect(getAnyMaster(scene)).toBe(root);

    const lateMaster = {
      name: 'bond_master_after',
      getClassName: () => 'Mesh',
      parent: root,
    };
    root.addChild(lateMaster);
    scene.meshes.push(lateMaster);

    const after = refreshMoleculeMasters(scene, { force: true });
    expect(after).toEqual([root]);
    expect(getMoleculeMasters(scene)).toEqual([root]);
  });
});


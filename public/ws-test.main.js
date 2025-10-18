// With Option A, stubs are loaded via <script> and exposed on window.proto
const pb = (typeof window!=='undefined' && window.proto && window.proto.fairchem && window.proto.fairchem.session) || {};
console.log('Vec3 available?', !!pb.Vec3);
// your page bootstrap here…
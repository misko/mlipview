import { create, toBinary } from '@bufbuild/protobuf';
import {
    ClientActionSchema,
    ClientAction_UserInteractionSchema,
    ClientAction_Start_SimType,
} from '../public/proto/fairchem_local_server2/session_pb.js';

function logCase(label, msg) {
    const bytes = toBinary(ClientActionSchema, msg);
    console.log(label, 'len', bytes.length);
    console.log(JSON.stringify(msg));
}

// Case A: nested created message (snake_case fields)
const uiA = create(ClientAction_UserInteractionSchema, {
    atomic_numbers: [8, 1, 1],
    positions: [0, 0, 0, 0, 0.757, 0.586, 0, -0.757, 0.586],
});
const msgA = create(ClientActionSchema, { seq: 1, schema_version: 1, user_interaction: uiA });
logCase('A', msgA);

// Case B: plain nested object (snake_case)
const msgB = create(ClientActionSchema, {
    seq: 2, schema_version: 1, user_interaction: {
        atomic_numbers: [8, 1, 1], positions: [0, 0, 0, 0, 0.757, 0.586, 0, -0.757, 0.586],
    }
});
logCase('B', msgB);

// Case C: start message (snake_case)
const msgC = create(ClientActionSchema, {
    seq: 3, schema_version: 1, start: {
        simulation_type: ClientAction_Start_SimType.MD,
    }
});
logCase('C', msgC);

// Case D: set via payload object with case/value
const msgD = create(ClientActionSchema, {
    seq: 4,
    schemaVersion: 1,
    payload: { case: 'userInteraction', value: { atomicNumbers: [8, 1, 1], positions: [0, 0, 0, 1, 0, 0, 0, 1, 0] } },
});
logCase('D', msgD);

// Case E: set via payload object with nested field name
const msgE = create(ClientActionSchema, {
    seq: 5,
    schemaVersion: 1,
    payload: { userInteraction: { atomicNumbers: [8, 1, 1], positions: [0, 0, 0, 1, 0, 0, 0, 1, 0] } },
});
logCase('E', msgE);

#include "common_dbc.h"

namespace {

const Signal sigs_880[] = {
    {
      .name = "EPAS_steeringRackForce",
      .b1 = 6,
      .b2 = 10,
      .bo = 48,
      .is_signed = false,
      .factor = 50,
      .offset = -25575.0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_steeringFault",
      .b1 = 5,
      .b2 = 1,
      .bo = 58,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_steeringReduced",
      .b1 = 4,
      .b2 = 1,
      .bo = 59,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_currentTuneMode",
      .b1 = 0,
      .b2 = 4,
      .bo = 60,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_torsionBarTorque",
      .b1 = 20,
      .b2 = 12,
      .bo = 32,
      .is_signed = false,
      .factor = 0.01,
      .offset = -20.5,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_eacErrorCode",
      .b1 = 16,
      .b2 = 4,
      .bo = 44,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_internalSAS",
      .b1 = 34,
      .b2 = 14,
      .bo = 16,
      .is_signed = false,
      .factor = 0.1,
      .offset = -819.200012,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_handsOnLevel",
      .b1 = 32,
      .b2 = 2,
      .bo = 30,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_sysStatusCounter",
      .b1 = 52,
      .b2 = 4,
      .bo = 8,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_eacStatus",
      .b1 = 48,
      .b2 = 3,
      .bo = 13,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
    {
      .name = "EPAS_sysStatusChecksum",
      .b1 = 56,
      .b2 = 8,
      .bo = 0,
      .is_signed = false,
      .factor = 1,
      .offset = 0,
      .is_little_endian = false,
      .type = SignalType::DEFAULT,
    },
};

const Msg msgs[] = {
  {
    .name = "EPAS_sysStatus",
    .address = 0x370,
    .size = 8,
    .num_sigs = ARRAYSIZE(sigs_880),
    .sigs = sigs_880,
  },
};

const Val vals[] = {
    {
      .name = "EPAS_currentTuneMode",
      .address = 0x370,
      .def_val = "1 DM_COMFORT 3 DM_SPORT 2 DM_STANDARD 0 FAIL_SAFE_DEFAULT 4 RWD_COMFORT 6 RWD_SPORT 5 RWD_STANDARD 7 UNAVAILABLE",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_eacErrorCode",
      .address = 0x370,
      .def_val = "14 EAC_ERROR_EPB_INHIBIT 3 EAC_ERROR_HANDS_ON 7 EAC_ERROR_HIGH_ANGLE_RATE_REQ 9 EAC_ERROR_HIGH_ANGLE_RATE_SAFETY 6 EAC_ERROR_HIGH_ANGLE_REQ 8 EAC_ERROR_HIGH_ANGLE_SAFETY 10 EAC_ERROR_HIGH_MMOT_SAFETY 11 EAC_ERROR_HIGH_TORSION_SAFETY 0 EAC_ERROR_IDLE 12 EAC_ERROR_LOW_ASSIST 2 EAC_ERROR_MAX_SPEED 1 EAC_ERROR_MIN_SPEED 13 EAC_ERROR_PINION_VEL_DIFF 4 EAC_ERROR_TMP_FAULT 5 EAR_ERROR_MAX_STEER_DELTA 15 SNA",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_eacStatus",
      .address = 0x370,
      .def_val = "2 EAC_ACTIVE 1 EAC_AVAILABLE 3 EAC_FAULT 0 EAC_INHIBITED 4 SNA",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_handsOnLevel",
      .address = 0x370,
      .def_val = "0 0 1 1 2 2 3 3",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_steeringFault",
      .address = 0x370,
      .def_val = "1 FAULT 0 NO_FAULT",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_steeringRackForce",
      .address = 0x370,
      .def_val = "1022 NOT_IN_SPEC 1023 SNA",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_steeringReduced",
      .address = 0x370,
      .def_val = "0 NORMAL_ASSIST 1 REDUCED_ASSIST",
      .sigs = sigs_880,
    },
    {
      .name = "EPAS_torsionBarTorque",
      .address = 0x370,
      .def_val = "0 SEE_SPECIFICATION 4095 SNA 4094 UNDEFINABLE_DATA",
      .sigs = sigs_880,
    },
};

}

const DBC tesla_can_epas = {
  .name = "tesla_can_epas",
  .num_msgs = ARRAYSIZE(msgs),
  .msgs = msgs,
  .vals = vals,
  .num_vals = ARRAYSIZE(vals),
};

dbc_init(tesla_can_epas)
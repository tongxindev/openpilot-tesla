#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2641264040911024235) {
   out_2641264040911024235[0] = delta_x[0] + nom_x[0];
   out_2641264040911024235[1] = delta_x[1] + nom_x[1];
   out_2641264040911024235[2] = delta_x[2] + nom_x[2];
   out_2641264040911024235[3] = delta_x[3] + nom_x[3];
   out_2641264040911024235[4] = delta_x[4] + nom_x[4];
   out_2641264040911024235[5] = delta_x[5] + nom_x[5];
   out_2641264040911024235[6] = delta_x[6] + nom_x[6];
   out_2641264040911024235[7] = delta_x[7] + nom_x[7];
   out_2641264040911024235[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8682440873081575844) {
   out_8682440873081575844[0] = -nom_x[0] + true_x[0];
   out_8682440873081575844[1] = -nom_x[1] + true_x[1];
   out_8682440873081575844[2] = -nom_x[2] + true_x[2];
   out_8682440873081575844[3] = -nom_x[3] + true_x[3];
   out_8682440873081575844[4] = -nom_x[4] + true_x[4];
   out_8682440873081575844[5] = -nom_x[5] + true_x[5];
   out_8682440873081575844[6] = -nom_x[6] + true_x[6];
   out_8682440873081575844[7] = -nom_x[7] + true_x[7];
   out_8682440873081575844[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6051135439158672014) {
   out_6051135439158672014[0] = 1.0;
   out_6051135439158672014[1] = 0;
   out_6051135439158672014[2] = 0;
   out_6051135439158672014[3] = 0;
   out_6051135439158672014[4] = 0;
   out_6051135439158672014[5] = 0;
   out_6051135439158672014[6] = 0;
   out_6051135439158672014[7] = 0;
   out_6051135439158672014[8] = 0;
   out_6051135439158672014[9] = 0;
   out_6051135439158672014[10] = 1.0;
   out_6051135439158672014[11] = 0;
   out_6051135439158672014[12] = 0;
   out_6051135439158672014[13] = 0;
   out_6051135439158672014[14] = 0;
   out_6051135439158672014[15] = 0;
   out_6051135439158672014[16] = 0;
   out_6051135439158672014[17] = 0;
   out_6051135439158672014[18] = 0;
   out_6051135439158672014[19] = 0;
   out_6051135439158672014[20] = 1.0;
   out_6051135439158672014[21] = 0;
   out_6051135439158672014[22] = 0;
   out_6051135439158672014[23] = 0;
   out_6051135439158672014[24] = 0;
   out_6051135439158672014[25] = 0;
   out_6051135439158672014[26] = 0;
   out_6051135439158672014[27] = 0;
   out_6051135439158672014[28] = 0;
   out_6051135439158672014[29] = 0;
   out_6051135439158672014[30] = 1.0;
   out_6051135439158672014[31] = 0;
   out_6051135439158672014[32] = 0;
   out_6051135439158672014[33] = 0;
   out_6051135439158672014[34] = 0;
   out_6051135439158672014[35] = 0;
   out_6051135439158672014[36] = 0;
   out_6051135439158672014[37] = 0;
   out_6051135439158672014[38] = 0;
   out_6051135439158672014[39] = 0;
   out_6051135439158672014[40] = 1.0;
   out_6051135439158672014[41] = 0;
   out_6051135439158672014[42] = 0;
   out_6051135439158672014[43] = 0;
   out_6051135439158672014[44] = 0;
   out_6051135439158672014[45] = 0;
   out_6051135439158672014[46] = 0;
   out_6051135439158672014[47] = 0;
   out_6051135439158672014[48] = 0;
   out_6051135439158672014[49] = 0;
   out_6051135439158672014[50] = 1.0;
   out_6051135439158672014[51] = 0;
   out_6051135439158672014[52] = 0;
   out_6051135439158672014[53] = 0;
   out_6051135439158672014[54] = 0;
   out_6051135439158672014[55] = 0;
   out_6051135439158672014[56] = 0;
   out_6051135439158672014[57] = 0;
   out_6051135439158672014[58] = 0;
   out_6051135439158672014[59] = 0;
   out_6051135439158672014[60] = 1.0;
   out_6051135439158672014[61] = 0;
   out_6051135439158672014[62] = 0;
   out_6051135439158672014[63] = 0;
   out_6051135439158672014[64] = 0;
   out_6051135439158672014[65] = 0;
   out_6051135439158672014[66] = 0;
   out_6051135439158672014[67] = 0;
   out_6051135439158672014[68] = 0;
   out_6051135439158672014[69] = 0;
   out_6051135439158672014[70] = 1.0;
   out_6051135439158672014[71] = 0;
   out_6051135439158672014[72] = 0;
   out_6051135439158672014[73] = 0;
   out_6051135439158672014[74] = 0;
   out_6051135439158672014[75] = 0;
   out_6051135439158672014[76] = 0;
   out_6051135439158672014[77] = 0;
   out_6051135439158672014[78] = 0;
   out_6051135439158672014[79] = 0;
   out_6051135439158672014[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8561185641601811804) {
   out_8561185641601811804[0] = state[0];
   out_8561185641601811804[1] = state[1];
   out_8561185641601811804[2] = state[2];
   out_8561185641601811804[3] = state[3];
   out_8561185641601811804[4] = state[4];
   out_8561185641601811804[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8561185641601811804[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8561185641601811804[7] = state[7];
   out_8561185641601811804[8] = state[8];
}
void F_fun(double *state, double dt, double *out_993737228106515281) {
   out_993737228106515281[0] = 1;
   out_993737228106515281[1] = 0;
   out_993737228106515281[2] = 0;
   out_993737228106515281[3] = 0;
   out_993737228106515281[4] = 0;
   out_993737228106515281[5] = 0;
   out_993737228106515281[6] = 0;
   out_993737228106515281[7] = 0;
   out_993737228106515281[8] = 0;
   out_993737228106515281[9] = 0;
   out_993737228106515281[10] = 1;
   out_993737228106515281[11] = 0;
   out_993737228106515281[12] = 0;
   out_993737228106515281[13] = 0;
   out_993737228106515281[14] = 0;
   out_993737228106515281[15] = 0;
   out_993737228106515281[16] = 0;
   out_993737228106515281[17] = 0;
   out_993737228106515281[18] = 0;
   out_993737228106515281[19] = 0;
   out_993737228106515281[20] = 1;
   out_993737228106515281[21] = 0;
   out_993737228106515281[22] = 0;
   out_993737228106515281[23] = 0;
   out_993737228106515281[24] = 0;
   out_993737228106515281[25] = 0;
   out_993737228106515281[26] = 0;
   out_993737228106515281[27] = 0;
   out_993737228106515281[28] = 0;
   out_993737228106515281[29] = 0;
   out_993737228106515281[30] = 1;
   out_993737228106515281[31] = 0;
   out_993737228106515281[32] = 0;
   out_993737228106515281[33] = 0;
   out_993737228106515281[34] = 0;
   out_993737228106515281[35] = 0;
   out_993737228106515281[36] = 0;
   out_993737228106515281[37] = 0;
   out_993737228106515281[38] = 0;
   out_993737228106515281[39] = 0;
   out_993737228106515281[40] = 1;
   out_993737228106515281[41] = 0;
   out_993737228106515281[42] = 0;
   out_993737228106515281[43] = 0;
   out_993737228106515281[44] = 0;
   out_993737228106515281[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_993737228106515281[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_993737228106515281[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_993737228106515281[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_993737228106515281[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_993737228106515281[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_993737228106515281[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_993737228106515281[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_993737228106515281[53] = -9.8000000000000007*dt;
   out_993737228106515281[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_993737228106515281[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_993737228106515281[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_993737228106515281[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_993737228106515281[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_993737228106515281[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_993737228106515281[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_993737228106515281[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_993737228106515281[62] = 0;
   out_993737228106515281[63] = 0;
   out_993737228106515281[64] = 0;
   out_993737228106515281[65] = 0;
   out_993737228106515281[66] = 0;
   out_993737228106515281[67] = 0;
   out_993737228106515281[68] = 0;
   out_993737228106515281[69] = 0;
   out_993737228106515281[70] = 1;
   out_993737228106515281[71] = 0;
   out_993737228106515281[72] = 0;
   out_993737228106515281[73] = 0;
   out_993737228106515281[74] = 0;
   out_993737228106515281[75] = 0;
   out_993737228106515281[76] = 0;
   out_993737228106515281[77] = 0;
   out_993737228106515281[78] = 0;
   out_993737228106515281[79] = 0;
   out_993737228106515281[80] = 1;
}
void h_25(double *state, double *unused, double *out_4562119136135831163) {
   out_4562119136135831163[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3607540865635187869) {
   out_3607540865635187869[0] = 0;
   out_3607540865635187869[1] = 0;
   out_3607540865635187869[2] = 0;
   out_3607540865635187869[3] = 0;
   out_3607540865635187869[4] = 0;
   out_3607540865635187869[5] = 0;
   out_3607540865635187869[6] = 1;
   out_3607540865635187869[7] = 0;
   out_3607540865635187869[8] = 0;
}
void h_24(double *state, double *unused, double *out_3023611308669197401) {
   out_3023611308669197401[0] = state[4];
   out_3023611308669197401[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8354050192601696903) {
   out_8354050192601696903[0] = 0;
   out_8354050192601696903[1] = 0;
   out_8354050192601696903[2] = 0;
   out_8354050192601696903[3] = 0;
   out_8354050192601696903[4] = 1;
   out_8354050192601696903[5] = 0;
   out_8354050192601696903[6] = 0;
   out_8354050192601696903[7] = 0;
   out_8354050192601696903[8] = 0;
   out_8354050192601696903[9] = 0;
   out_8354050192601696903[10] = 0;
   out_8354050192601696903[11] = 0;
   out_8354050192601696903[12] = 0;
   out_8354050192601696903[13] = 0;
   out_8354050192601696903[14] = 1;
   out_8354050192601696903[15] = 0;
   out_8354050192601696903[16] = 0;
   out_8354050192601696903[17] = 0;
}
void h_30(double *state, double *unused, double *out_5387964256309980461) {
   out_5387964256309980461[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3309149475856428886) {
   out_3309149475856428886[0] = 0;
   out_3309149475856428886[1] = 0;
   out_3309149475856428886[2] = 0;
   out_3309149475856428886[3] = 0;
   out_3309149475856428886[4] = 1;
   out_3309149475856428886[5] = 0;
   out_3309149475856428886[6] = 0;
   out_3309149475856428886[7] = 0;
   out_3309149475856428886[8] = 0;
}
void h_26(double *state, double *unused, double *out_243074193194332924) {
   out_243074193194332924[0] = state[7];
}
void H_26(double *state, double *unused, double *out_303014895874387268) {
   out_303014895874387268[0] = 0;
   out_303014895874387268[1] = 0;
   out_303014895874387268[2] = 0;
   out_303014895874387268[3] = 0;
   out_303014895874387268[4] = 0;
   out_303014895874387268[5] = 0;
   out_303014895874387268[6] = 0;
   out_303014895874387268[7] = 1;
   out_303014895874387268[8] = 0;
}
void h_27(double *state, double *unused, double *out_7978392425784901179) {
   out_7978392425784901179[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1134386164056003975) {
   out_1134386164056003975[0] = 0;
   out_1134386164056003975[1] = 0;
   out_1134386164056003975[2] = 0;
   out_1134386164056003975[3] = 1;
   out_1134386164056003975[4] = 0;
   out_1134386164056003975[5] = 0;
   out_1134386164056003975[6] = 0;
   out_1134386164056003975[7] = 0;
   out_1134386164056003975[8] = 0;
}
void h_29(double *state, double *unused, double *out_3697516421574299501) {
   out_3697516421574299501[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3819380820170821070) {
   out_3819380820170821070[0] = 0;
   out_3819380820170821070[1] = 1;
   out_3819380820170821070[2] = 0;
   out_3819380820170821070[3] = 0;
   out_3819380820170821070[4] = 0;
   out_3819380820170821070[5] = 0;
   out_3819380820170821070[6] = 0;
   out_3819380820170821070[7] = 0;
   out_3819380820170821070[8] = 0;
}
void h_28(double *state, double *unused, double *out_2976719158124055534) {
   out_2976719158124055534[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1263018196898709504) {
   out_1263018196898709504[0] = 1;
   out_1263018196898709504[1] = 0;
   out_1263018196898709504[2] = 0;
   out_1263018196898709504[3] = 0;
   out_1263018196898709504[4] = 0;
   out_1263018196898709504[5] = 0;
   out_1263018196898709504[6] = 0;
   out_1263018196898709504[7] = 0;
   out_1263018196898709504[8] = 0;
}
void h_31(double *state, double *unused, double *out_6003352286819671268) {
   out_6003352286819671268[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3469134384876629384) {
   out_3469134384876629384[0] = 0;
   out_3469134384876629384[1] = 0;
   out_3469134384876629384[2] = 0;
   out_3469134384876629384[3] = 0;
   out_3469134384876629384[4] = 0;
   out_3469134384876629384[5] = 0;
   out_3469134384876629384[6] = 0;
   out_3469134384876629384[7] = 0;
   out_3469134384876629384[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_2641264040911024235) {
  err_fun(nom_x, delta_x, out_2641264040911024235);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8682440873081575844) {
  inv_err_fun(nom_x, true_x, out_8682440873081575844);
}
void car_H_mod_fun(double *state, double *out_6051135439158672014) {
  H_mod_fun(state, out_6051135439158672014);
}
void car_f_fun(double *state, double dt, double *out_8561185641601811804) {
  f_fun(state,  dt, out_8561185641601811804);
}
void car_F_fun(double *state, double dt, double *out_993737228106515281) {
  F_fun(state,  dt, out_993737228106515281);
}
void car_h_25(double *state, double *unused, double *out_4562119136135831163) {
  h_25(state, unused, out_4562119136135831163);
}
void car_H_25(double *state, double *unused, double *out_3607540865635187869) {
  H_25(state, unused, out_3607540865635187869);
}
void car_h_24(double *state, double *unused, double *out_3023611308669197401) {
  h_24(state, unused, out_3023611308669197401);
}
void car_H_24(double *state, double *unused, double *out_8354050192601696903) {
  H_24(state, unused, out_8354050192601696903);
}
void car_h_30(double *state, double *unused, double *out_5387964256309980461) {
  h_30(state, unused, out_5387964256309980461);
}
void car_H_30(double *state, double *unused, double *out_3309149475856428886) {
  H_30(state, unused, out_3309149475856428886);
}
void car_h_26(double *state, double *unused, double *out_243074193194332924) {
  h_26(state, unused, out_243074193194332924);
}
void car_H_26(double *state, double *unused, double *out_303014895874387268) {
  H_26(state, unused, out_303014895874387268);
}
void car_h_27(double *state, double *unused, double *out_7978392425784901179) {
  h_27(state, unused, out_7978392425784901179);
}
void car_H_27(double *state, double *unused, double *out_1134386164056003975) {
  H_27(state, unused, out_1134386164056003975);
}
void car_h_29(double *state, double *unused, double *out_3697516421574299501) {
  h_29(state, unused, out_3697516421574299501);
}
void car_H_29(double *state, double *unused, double *out_3819380820170821070) {
  H_29(state, unused, out_3819380820170821070);
}
void car_h_28(double *state, double *unused, double *out_2976719158124055534) {
  h_28(state, unused, out_2976719158124055534);
}
void car_H_28(double *state, double *unused, double *out_1263018196898709504) {
  H_28(state, unused, out_1263018196898709504);
}
void car_h_31(double *state, double *unused, double *out_6003352286819671268) {
  h_31(state, unused, out_6003352286819671268);
}
void car_H_31(double *state, double *unused, double *out_3469134384876629384) {
  H_31(state, unused, out_3469134384876629384);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);

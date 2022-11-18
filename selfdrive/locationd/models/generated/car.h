#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_2641264040911024235);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8682440873081575844);
void car_H_mod_fun(double *state, double *out_6051135439158672014);
void car_f_fun(double *state, double dt, double *out_8561185641601811804);
void car_F_fun(double *state, double dt, double *out_993737228106515281);
void car_h_25(double *state, double *unused, double *out_4562119136135831163);
void car_H_25(double *state, double *unused, double *out_3607540865635187869);
void car_h_24(double *state, double *unused, double *out_3023611308669197401);
void car_H_24(double *state, double *unused, double *out_8354050192601696903);
void car_h_30(double *state, double *unused, double *out_5387964256309980461);
void car_H_30(double *state, double *unused, double *out_3309149475856428886);
void car_h_26(double *state, double *unused, double *out_243074193194332924);
void car_H_26(double *state, double *unused, double *out_303014895874387268);
void car_h_27(double *state, double *unused, double *out_7978392425784901179);
void car_H_27(double *state, double *unused, double *out_1134386164056003975);
void car_h_29(double *state, double *unused, double *out_3697516421574299501);
void car_H_29(double *state, double *unused, double *out_3819380820170821070);
void car_h_28(double *state, double *unused, double *out_2976719158124055534);
void car_H_28(double *state, double *unused, double *out_1263018196898709504);
void car_h_31(double *state, double *unused, double *out_6003352286819671268);
void car_H_31(double *state, double *unused, double *out_3469134384876629384);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}
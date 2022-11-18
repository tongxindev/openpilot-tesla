#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_3155115096055378676);
void live_err_fun(double *nom_x, double *delta_x, double *out_8395496394702832161);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8916956029375892266);
void live_H_mod_fun(double *state, double *out_4679673059438427241);
void live_f_fun(double *state, double dt, double *out_8106099418900950754);
void live_F_fun(double *state, double dt, double *out_4673684726810862167);
void live_h_4(double *state, double *unused, double *out_5235680874440460819);
void live_H_4(double *state, double *unused, double *out_3623290639705711408);
void live_h_9(double *state, double *unused, double *out_4479381555540235079);
void live_H_9(double *state, double *unused, double *out_3864480286335302053);
void live_h_10(double *state, double *unused, double *out_2164176216250742354);
void live_H_10(double *state, double *unused, double *out_3212870814963671173);
void live_h_12(double *state, double *unused, double *out_4049753783849693471);
void live_H_12(double *state, double *unused, double *out_8642747047737673203);
void live_h_31(double *state, double *unused, double *out_8663581636506766747);
void live_H_31(double *state, double *unused, double *out_7058433993646864704);
void live_h_32(double *state, double *unused, double *out_1427867763541487384);
void live_H_32(double *state, double *unused, double *out_4274775829228447111);
void live_h_13(double *state, double *unused, double *out_2549469855511810888);
void live_H_13(double *state, double *unused, double *out_7089291747604569330);
void live_h_14(double *state, double *unused, double *out_4479381555540235079);
void live_H_14(double *state, double *unused, double *out_3864480286335302053);
void live_h_33(double *state, double *unused, double *out_2749875488518733959);
void live_H_33(double *state, double *unused, double *out_3907876989008007100);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}
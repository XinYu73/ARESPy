/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "config.h"

#ifdef HAVE_FORTRAN

#include "xc.h"
#include "string_f.h"

/* xc_config.h needs to be included to use FLOAT and related macros*/
#include "xc_config.h"

/* version */
void XC_FC_FUNC(f90_version, F90_VERSION)
     (int *major, int *minor, int *micro)
{
  XC(version)(major, minor, micro);
}

void XC_FC_FUNC(f90_version_string, F90_VERSIN_STRING)
     (STR_F_TYPE version_string STR_ARG1)
{
  const char *version;

  version = XC(version_string)();
  TO_F_STR1(version, version_string);
}

/* info */
CC_FORTRAN_INT XC_FC_FUNC(f90_info_number, F90_INFO_NUMBER)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->number;
}


CC_FORTRAN_INT XC_FC_FUNC(f90_info_kind, F90_INFO_KIND)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->kind;
}


void XC_FC_FUNC(f90_info_name, F90_INFO_NAME)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((XC(func_info_type) *)(*info))->name, s);
}


CC_FORTRAN_INT  XC_FC_FUNC(f90_info_family, F90_INFO_FAMILY)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->family;
}


CC_FORTRAN_INT  XC_FC_FUNC(f90_info_flags, F90_INFO_FLAGS)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->flags;
}


void XC_FC_FUNC(f90_info_refs, F90_INFO_REFS)
     (void **info, CC_FORTRAN_INT *number, STR_F_TYPE ref_f STR_ARG1)
{
  XC(func_info_type) *func_p = (XC(func_info_type) *)(*info);

  assert(*number >=0 && *number < 5);

  if(func_p->refs[*number] == NULL){
    *number = -1;
    return;
  }

  TO_F_STR1(func_p->refs[*number]->ref, ref_f);

  (*number)++;
  fflush(stdout);
}


void XC_FC_FUNC(f90_functional_get_name, F90_FUNCTIONAL_GET_NAME)
     (CC_FORTRAN_INT *func_number, STR_F_TYPE func_string STR_ARG1)
{
  char *name;

  name = XC(functional_get_name)(*func_number);
  if ( name == NULL ) name = strdup("unknown");

  TO_F_STR1(name, func_string);
  free(name);
}


CC_FORTRAN_INT  XC_FC_FUNC(f90_functional_get_number, F90_FUNCTIONAL_GET_NUMBER)
     (STR_F_TYPE func_string STR_ARG1)
{
  char *name;
  int ret;

  TO_C_STR1(func_string, name);
  
  ret = XC(functional_get_number)(name);
  free(name);

  return (CC_FORTRAN_INT) ret;
}


/* functionals */
CC_FORTRAN_INT  XC_FC_FUNC(f90_family_from_id, F90_FAMILY_FROM_ID)
  (CC_FORTRAN_INT  *functional)
{
  return (CC_FORTRAN_INT) XC(family_from_id)((int) (*functional), NULL, NULL);
}


/* Standard initialization */
void XC_FC_FUNC(f90_func_init, F90_FUNC_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  XC(func_type) *func_p;
  
  func_p = (XC(func_type) *)malloc(sizeof(XC(func_type)));
  XC(func_init)(func_p, (int) (*functional), (int) (*nspin));

  *p    = (void *) func_p;
  *info = (void *)(func_p->info);
}


void XC_FC_FUNC(f90_func_end, F90_FUNC_END)
     (void **p)
{
  XC(func_end)((XC(func_type) *)(*p));
  free(*p);
  *p = NULL;
}


/* LDAs */

void XC_FC_FUNC(f90_lda, F90_LDA)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, 
      FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, zk, vrho, v2rho2, v3rho3);
}

void XC_FC_FUNC(f90_lda_exc, F90_LDA_EXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho,
      FLOAT *zk)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, zk, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_lda_exc_vxc, F90_LDA_EXC_VXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, 
      FLOAT *zk, FLOAT *vrho)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, zk, vrho, NULL, NULL);
}

void XC_FC_FUNC(f90_lda_vxc, F90_LDA_VXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, 
      FLOAT *vrho)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, NULL, vrho, NULL, NULL);
}

void XC_FC_FUNC(f90_lda_vxc_fxc, F90_LDA_VXC_FXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho,
      FLOAT *vrho, FLOAT *v2rho2)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, NULL, vrho, v2rho2, NULL);
}

void XC_FC_FUNC(f90_lda_fxc, F90_LDA_FXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho,
      FLOAT *v2rho2)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, NULL, NULL, v2rho2, NULL);
}

void XC_FC_FUNC(f90_lda_kxc, F90_LDA_KXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho,
      FLOAT *v3rho3)
{
  XC(lda)((XC(func_type) *)(*p), *np, rho, NULL, NULL, NULL, v3rho3);
}


/* GGAs */

void XC_FC_FUNC(f90_gga, F90_GGA)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2,
      FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, zk, vrho, vsigma, 
	  v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);
}

void XC_FC_FUNC(f90_gga_exc, F90_GGA_EXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, zk, NULL, NULL, 
	  NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_gga_exc_vxc, F90_GGA_EXC_VXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, zk, vrho, vsigma, 
	  NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_gga_vxc, F90_GGA_VXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, NULL, vrho, vsigma, 
	  NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_gga_vxc_fxc, F90_GGA_VXC_FXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *vrho, FLOAT *vsigma,
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, NULL, vrho, vsigma, 
	  v2rho2, v2rhosigma, v2sigma2, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_gga_fxc, F90_GGA_FXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, NULL, NULL, NULL, 
	  v2rho2, v2rhosigma, v2sigma2, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_gga_kxc, F90_GGA_KXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, 
      FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3)
{
  XC(gga)((XC(func_type) *)(*p), *np, rho, sigma, NULL, NULL, NULL, 
	  NULL, NULL, NULL, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);
}

void XC_FC_FUNC(f90_nlc_coef, F90_nlc_COEF)
  (void **p, FLOAT *nlc_b, FLOAT *nlc_c)
{
  XC(nlc_coef)((XC(func_type) *)(*p), nlc_b, nlc_c);
}

void XC_FC_FUNC(f90_gga_lb_modified, F90_GGA_LB_MODIFIED)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *r, FLOAT *vrho)
{
  XC(gga_lb_modified)((XC(func_type) *)(*p), *np, rho, sigma, *r, vrho);
}

void XC_FC_FUNC(f90_gga_ak13_get_asymptotic, F90_GGA_AK13_GET_ASYMPTOTIC)
  (FLOAT *homo, FLOAT *asymp)
{
  *asymp = XC(gga_ak13_get_asymptotic)(*homo);
}

void XC_FC_FUNC(f90_hyb_exx_coef, F90_HYB_EXX_COEF)
   (void **p, FLOAT *coef)
{
  *coef = XC(hyb_exx_coef)((XC(func_type) *)(*p));
}

void XC_FC_FUNC(f90_hyb_cam_coef, F90_HYB_CAM_COEF)
  (void **p, FLOAT *omega, FLOAT *alpha, FLOAT *beta)
{
  XC(hyb_cam_coef)((XC(func_type) *)(*p), omega, alpha, beta);
}


/* meta-GGAs */

void XC_FC_FUNC(f90_mgga, F90_MGGA)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *lapl, FLOAT *tau,
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau,
      FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
      FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
      FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau)
{
  XC(mgga)((XC(func_type) *)(*p), *np, rho, sigma, lapl, tau, 
	   zk, vrho, vsigma, vlapl, vtau,
	   v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, 
	   v2sigmalapl, v2sigmatau, v2lapltau);

}

void XC_FC_FUNC(f90_mgga_exc, F90_MGGA_EXC)
     (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *lapl, FLOAT *tau, 
      FLOAT *zk)
{
  XC(mgga)((XC(func_type) *)(*p), *np, rho, sigma, lapl, tau, 
	   zk, NULL, NULL, NULL, NULL, 
	   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_mgga_exc_vxc, F90_MGGA_EXC_VXC)
  (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *lapl, FLOAT *tau,
   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau)
{
  XC(mgga)((XC(func_type) *)(*p), *np, rho, sigma, lapl, tau, 
	   zk, vrho, vsigma, vlapl, vtau, 
	   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_mgga_vxc, F90_MGGA_VXC)
  (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *lapl, FLOAT *tau,
   FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau)
{
  XC(mgga)((XC(func_type) *)(*p), *np, rho, sigma, lapl, tau, 
	   NULL, vrho, vsigma, vlapl, vtau, 
	   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_mgga_vxc_fxc, F90_MGGA_VXC_FXC)
  (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *lapl, FLOAT *tau,
   FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau,
   FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
   FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
   FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau)
{
  XC(mgga)((XC(func_type) *)(*p), *np, rho, sigma, lapl, tau, 
	   NULL, vrho, vsigma, vlapl, vtau,
	   v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, 
	   v2sigmalapl, v2sigmatau, v2lapltau);
}

void XC_FC_FUNC(f90_mgga_fxc, F90_MGGA_FXC)
  (void **p, CC_FORTRAN_INT *np, FLOAT *rho, FLOAT *sigma, FLOAT *lapl, FLOAT *tau,
      FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
      FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
      FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau)
{
  XC(mgga)((XC(func_type) *)(*p), *np, rho, sigma, lapl, tau, 
	   NULL, NULL, NULL, NULL, NULL, 
	   v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, 
	   v2sigmalapl, v2sigmatau, v2lapltau);
}


void XC_FC_FUNC(f90_lda_c_1d_csc_set_par, F90_LDA_C_1D_CSC_SET_PAR)
  (void **p, CC_FORTRAN_INT *interaction, FLOAT *bb)
{
  XC(lda_c_1d_csc_set_params)((XC(func_type) *)(*p), *interaction, *bb);
}

void XC_FC_FUNC(f90_lda_c_2d_prm_set_par, F90_LDA_C_2D_PRM_SET_PAR)
  (void **p, FLOAT *N)
{
  XC(lda_c_2d_prm_set_params)((XC(func_type) *)(*p), *N);
}

void XC_FC_FUNC(f90_lda_c_vwn_set_par, F90_LDA_C_VWN_SET_PAR)
  (void **p, CC_FORTRAN_INT *spin_interpolation)
{
  XC(lda_c_vwn_set_params)((XC(func_type) *)(*p), *spin_interpolation);
}

void XC_FC_FUNC(f90_lda_c_xalpha_set_par, F90_LDA_C_XALPHA_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(lda_c_xalpha_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_lda_x_set_par, F90_LDA_X_SET_PAR)
  (void **p, FLOAT *alpha, CC_FORTRAN_INT *relativistic, FLOAT *omega)
{
  XC(lda_x_set_params)((XC(func_type) *)(*p), *alpha, *relativistic, *omega);
}

void XC_FC_FUNC(f90_lda_x_1d_set_par, F90_LDA_X_1D_SET_PAR)
  (void **p, CC_FORTRAN_INT *interaction, FLOAT *bb)
{
  XC(lda_x_1d_set_params)((XC(func_type) *)(*p), *interaction, *bb);
}

void XC_FC_FUNC(f90_lda_xc_ksdt_set_par, F90_LDA_XC_KSDT_SET_PAR)
  (void **p, FLOAT *T)
{
  XC(lda_xc_ksdt_set_params)((XC(func_type) *)(*p), *T);
}

void XC_FC_FUNC(f90_gga_c_lyp_set_par, F90_GGA_C_LYP_SET_PAR)
  (void **p, FLOAT *A, FLOAT *B, FLOAT *c, FLOAT *d)
{
  XC(gga_c_lyp_set_params)((XC(func_type) *)(*p), *A, *B, *c, *d);
}

void XC_FC_FUNC(f90_gga_c_pbe_set_par, F90_GGA_C_PBE_SET_PAR)
  (void **p, FLOAT *beta)
{
  XC(gga_c_pbe_set_params)((XC(func_type) *)(*p), *beta);
}

void XC_FC_FUNC(f90_gga_k_tflw_set_par, F90_GGA_K_TFLW_SET_PAR)
  (void **p, FLOAT *gamma, FLOAT *lambda, FLOAT *N)
{
  XC(gga_k_tflw_set_params)((XC(func_type) *)(*p), *gamma, *lambda, *N);
}

void XC_FC_FUNC(f90_gga_x_2d_b88_set_par, F90_GGA_X_2D_B88_SET_PAR)
  (void **p, FLOAT *beta)
{
  XC(gga_x_2d_b88_set_params)((XC(func_type) *)(*p), *beta);
}

void XC_FC_FUNC(f90_gga_x_b86_set_par, F90_GGA_X_B86_SET_PAR)
  (void **p, FLOAT *beta, FLOAT *gamma, FLOAT *omega)
{
  XC(gga_x_b86_set_params)((XC(func_type) *)(*p), *beta, *gamma, *omega);
}

void XC_FC_FUNC(f90_gga_x_b88_set_par, F90_GGA_X_B88_SET_PAR)
  (void **p, FLOAT *beta, FLOAT *gamma)
{
  XC(gga_x_b88_set_params)((XC(func_type) *)(*p), *beta, *gamma);
}

void XC_FC_FUNC(f90_gga_x_hjs_set_par, F90_GGA_X_HJS_SET_PAR)
  (void **p, FLOAT *omega)
{
  XC(gga_x_hjs_set_params)((XC(func_type) *)(*p), *omega);
}

void XC_FC_FUNC(f90_gga_x_ityh_set_par, F90_GGA_X_ITYH_SET_PAR)
  (void **p, CC_FORTRAN_INT *func_id, FLOAT *omega)
{
  XC(gga_x_ityh_set_params)((XC(func_type) *)(*p), *func_id, *omega);
}

void XC_FC_FUNC(f90_gga_x_kt_set_par, F90_GGA_X_KT_SET_PAR)
  (void **p, FLOAT *gamma, FLOAT *delta)
{
  XC(gga_x_kt_set_params)((XC(func_type) *)(*p), *gamma, *delta);
}

void XC_FC_FUNC(f90_gga_lb_set_par, F90_GGA_LB_SET_PAR)
  (void **p, CC_FORTRAN_INT *modified, FLOAT *threshold, FLOAT *ip, FLOAT *qtot)
{
  XC(gga_lb_set_params)((XC(func_type) *)(*p), *modified, *threshold, *ip, *qtot);
}

void XC_FC_FUNC(f90_gga_x_optx_set_par, F90_GGA_X_OPTX_SET_PAR)
  (void **p, FLOAT *a, FLOAT *b, FLOAT *gamma)
{
  XC(gga_x_optx_set_params)((XC(func_type) *)(*p), *a, *b, *gamma);
}

void XC_FC_FUNC(f90_gga_x_pbe_set_par, F90_GGA_X_PBE_SET_PAR)
  (void **p, FLOAT *kappa, FLOAT *mu)
{
  XC(gga_x_pbe_set_params)((XC(func_type) *)(*p), *kappa, *mu);
}

void XC_FC_FUNC(f90_gga_x_lambda_set_par, F90_GGA_X_LAMBDA_SET_PAR)
  (void **p, FLOAT *N)
{
  XC(gga_x_lambda_set_params)((XC(func_type) *)(*p), *N);
}

void XC_FC_FUNC(f90_gga_x_pw91_set_par, F90_GGA_X_PW91_SET_PAR)
  (void **p, FLOAT *a, FLOAT *b, FLOAT *c, FLOAT *d, FLOAT *f, FLOAT *alpha, FLOAT *expo)
{
  XC(gga_x_pw91_set_params)((XC(func_type) *)(*p), *a, *b, *c, *d, *f, *alpha, *expo);
}

void XC_FC_FUNC(f90_gga_x_pw91_set_par2, F90_GGA_X_PW91_SET_PAR2)
  (void **p, FLOAT *bt, FLOAT *alpha, FLOAT *expo)
{
  XC(gga_x_pw91_set_params2)((XC(func_type) *)(*p), *bt, *alpha, *expo);
}

void XC_FC_FUNC(f90_gga_x_rpbe_set_par, F90_GGA_X_RPBE_SET_PAR)
  (void **p, FLOAT *kappa, FLOAT *mu)
{
  XC(gga_x_rpbe_set_params)((XC(func_type) *)(*p), *kappa, *mu);
}

void XC_FC_FUNC(f90_gga_x_sfat_set_par, F90_GGA_X_SFAT_SET_PAR)
  (void **p, CC_FORTRAN_INT *func_id, FLOAT *omega)
{
  XC(gga_x_sfat_set_params)((XC(func_type) *)(*p), *func_id, *omega);
}

void XC_FC_FUNC(f90_gga_x_ssb_sw_set_par, F90_GGA_X_SSB_SW_SET_PAR)
  (void **p, FLOAT *A, FLOAT *B, FLOAT *C, FLOAT *D, FLOAT *E)
{
  XC(gga_x_ssb_sw_set_params)((XC(func_type) *)(*p), *A, *B, *C, *D, *E);
}

void XC_FC_FUNC(f90_gga_x_wpbeh_set_par, F90_GGA_X_WPBEH_SET_PAR)
  (void **p, FLOAT *omega)
{
  XC(gga_x_wpbeh_set_params)((XC(func_type) *)(*p), *omega);
}

void XC_FC_FUNC(f90_hyb_gga_xc_hse_set_par, F90_HYB_GGA_XC_HSE_SET_PAR)
  (void **p, FLOAT *beta, FLOAT *omega)
{
  XC(hyb_gga_xc_hse_set_params)((XC(func_type) *)(*p), *beta, *omega);
}

void XC_FC_FUNC(f90_hyb_gga_xc_pbeh_set_par, F90_HYB_GGA_XC_PBEH_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(hyb_gga_xc_pbeh_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_hyb_mgga_x_mvsh_set_par, F90_HYB_MGGA_X_MVSH_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(hyb_mgga_x_mvsh_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_hyb_mgga_xc_tpssh_set_par, F90_HYB_MGGA_XC_TPSSH_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(hyb_mgga_xc_tpssh_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_hyb_mgga_xc_revtpssh_set_par, F90_HYB_MGGA_XC_REVTPSSH_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(hyb_mgga_xc_revtpssh_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_mgga_c_bc95_set_par, F90_MGGA_C_BC95_SET_PAR)
  (void **p, FLOAT *css, FLOAT *copp)
{
  XC(mgga_c_bc95_set_params)((XC(func_type) *)(*p), *css, *copp);
}

void XC_FC_FUNC(f90_mgga_c_pkzb_set_par, F90_MGGA_C_PKZB_SET_PAR)
  (void **p, FLOAT *beta, FLOAT *d, FLOAT *C0_0, FLOAT *C0_1, FLOAT *C0_2, FLOAT *C0_3)
{
  XC(mgga_c_pkzb_set_params)((XC(func_type) *)(*p), *beta, *d, *C0_0, *C0_1, *C0_2, *C0_3);
}

void XC_FC_FUNC(f90_mgga_x_tb09_set_par, F90_MGGA_X_TB09_SET_PAR)
  (void **p, FLOAT *c)
{
  XC(mgga_x_tb09_set_params)((XC(func_type) *)(*p), *c);
}

void XC_FC_FUNC(f90_hyb_mgga_x_ms2h_set_par, F90_HYB_MGGA_X_MS2H_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(hyb_mgga_x_ms2h_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_hyb_mgga_x_scan0_set_par, F90_HYB_MGGA_X_SCAN0_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(hyb_mgga_x_scan0_set_params)((XC(func_type) *)(*p), *alpha);
}

void XC_FC_FUNC(f90_mgga_x_tpss_set_par, F90_MGGA_X_TPSS_SET_PAR)
  (void **p, FLOAT *b, FLOAT *c, FLOAT *e, FLOAT *kappa, FLOAT *mu)
{
  XC(mgga_x_tpss_set_params)((XC(func_type) *)(*p), *b, *c, *e, *kappa, *mu);
}

#endif

#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>


extern "C" {

    float directional_coupler_lc(float wavelength_nm, std::complex<float> n_eff_1, std::complex<float> n_eff_2, float& dn_eff, float& lc_m, float& lc_um);
    float grating_coupler_period(float wavelength, float& n_eff, float& n_clad, float& incidence_angle_deg, float& diffration_order);
    float loss(float wavelength, std::complex<float> n);
    float qpm_wavenumber(float pmp_n, float& pmp_l, float& sig_n, float&  sig_l, float&  idl_n, float&  idl_l, float& period_qpm, std::string type, int sgn_1, int sgn_2);
    float qpm_period(float pmp_n, float& pmp_l, float& sig_n, float&  sig_l, float&  idl_n, float&  idl_l, std::string type, int sgn_1, int sgn_2);

    double degToRad(double degrees) {
        double radians = degrees * M_PI / 180.0;
        return radians;
    }

    float directional_coupler_lc(float wavelength_nm, std::complex<float> n_eff_1, std::complex<float> n_eff_2, float& dn_eff, float& lc_m, float& lc_um) {
        float wavelength_m = wavelength_nm * 1.e-9;
        dn_eff = (n_eff_1 - n_eff_2).real();
        lc_m = wavelength_m / (2. * dn_eff);
        lc_um = lc_m * 1.e6; 
        return lc_um;
    }

    float grating_coupler_period(float wavelength, float& n_eff, float& n_clad, float& incidence_angle_deg, float& diffration_order){
        double pi = M_PI;
        float k0 = 2.0 * pi / wavelength;
        float beta = n_eff * k0;
        double radians = degToRad(incidence_angle_deg);
        float grating_period = (2.* pi *diffration_order) / (beta - k0*n_clad*sin(radians));
        return grating_period;
        // diffration_order=1 
    }

    float loss(float wavelength, std::complex<float> n){
        double pi = M_PI;
        float kappa = n.imag();
        float alpha = 4.34 * 4 * pi * abs(kappa) / wavelength;
        return alpha;
    }

    float qpm_wavenumber(float pmp_n, float& pmp_l, float& sig_n, float&  sig_l, float&  idl_n, float&  idl_l, float& period_qpm, std::string type, int sgn_1, int sgn_2){
        double pi = M_PI;
        float pi2 = pi * 2;
        float k_pmp = pmp_n * pi2 / pmp_l;
        float k_sig = sig_n * pi2 / sig_l;
        float k_idl = idl_n * pi2 / idl_l;
        float k_qpm = pi2 / period_qpm;
        if (type == "forward"){
            int sgn_1 = 1;
            int sgn_2 = 1;
        } else if (type == "forward_backward"){
            int sgn_1 = 1;
            int sgn_2 = -1;
        } else if (type == "backward"){
            int sgn_1 = -1;
            int sgn_2 = -1;
        }
        float k_mismatch = k_idl * sgn_1 + k_sig * sgn_2 + k_qpm - k_pmp;
        
        return k_mismatch;
    }

    float qpm_period(float pmp_n, float& pmp_l, float& sig_n, float&  sig_l, float&  idl_n, float&  idl_l, std::string type, int sgn_1, int sgn_2){

        double pi = M_PI;
        float pi2 = pi * 2;
        float k_pmp = pmp_n * pi2 / pmp_l;
        float k_sig = sig_n * pi2 / sig_l;
        float k_idl = idl_n * pi2 / idl_l;
        
        if (type == "forward"){
            int sgn_1 = 1;
            int sgn_2 = 1;
        } else if (type == "forward_backward") {
            int sgn_1 = 1;
            int sgn_2 = -1;
        } else if (type == "backward"){
            int sgn_1 = -1;
            int sgn_2 = -1;
        }
        
        float k_qpm = k_pmp - k_idl * sgn_1 - k_sig * sgn_2;
        float l_qpm = pi2 / k_qpm;
        return l_qpm;
    } 

}

int main() {
    return 0;
}
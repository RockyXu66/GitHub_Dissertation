function metric_par = hdrvdp_parse_options( options )
% HDRVDP_PARSE_OPTIONS (internal) parse HDR-VDP options and create two
% structures: view_cond with viewing conditions and metric_par with metric
% parameters
%
% Copyright (c) 2011, Rafal Mantiuk <mantiuk@gmail.com>

% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

% Defaults

% Peak contrast from Daly's CSF for L_adapt = 30 cd/m^2
daly_peak_contrast_sens = 0.006894596;

metric_par.sensitivity_correction = daly_peak_contrast_sens / 10.^-3.43304; 

metric_par.view_dist = 0.5;

metric_par.spectral_emission = [];


metric_par.orient_count = 4; % the number of orientations to consider

% Various optional features
metric_par.do_masking = true;
metric_par.do_mtf = true;
metric_par.do_spatial_pooling = true;
metric_par.noise_model = true;
metric_par.do_quality_raw_data = false; % for development purposes only


metric_par.steerpyr_filter = 'sp3Filters';

metric_par.mask_p = 0.544068;
metric_par.mask_self = 0.189065;
metric_par.mask_xo = 0.449199;
metric_par.mask_xn = 1.52512;
metric_par.mask_q = 0.49576;

metric_par.psych_func_slope = log10(3.5);
metric_par.beta = metric_par.psych_func_slope-metric_par.mask_p;

% Spatial summation
metric_par.si_slope = -0.850147;
metric_par.si_sigma = -0.000502005;
metric_par.si_ampl = 0;

% Cone and rod cvi functions
metric_par.cvi_sens_drop = 0.0704457;
metric_par.cvi_trans_slope = 0.0626528;
metric_par.cvi_low_slope = -0.00222585;

metric_par.rod_sensitivity = 0;
%metric_par.rod_sensitivity = -0.383324;
metric_par.cvi_sens_drop_rod = -0.58342;

% Achromatic CSF
metric_par.csf_m1_f_max = 0.425509;
metric_par.csf_m1_s_high = -0.227224;
metric_par.csf_m1_s_low = -0.227224;
metric_par.csf_m1_exp_low = log10( 2 );

% Daly CSF model
metric_par.csf_stim_area = 0; %2.5^2;
metric_par.csf_epsilon = -0.546385;
metric_par.csf_peak_shift = 0.235954;
metric_par.csf_lf_slope = -0.844601;
metric_par.csf_peak_shift_lum = 1.16336;
metric_par.csf_peak_shift_slope = -0.912733;

% Fix for the non-linearity after cortex transform
metric_par.ibf_fix = log10(0.321678);

% Rod CSF
metric_par.csf_rod_f_max = 0.15;
metric_par.csf_rod_s_low = -0.266858;    
metric_par.csf_rod_exp_low = log10(2);    
metric_par.csf_rod_s_high = -0.266858;    

metric_par.csf_params = [ ...
    0.3631    0.8320    2.8475    1.4717    4.3220
    1.0012    0.4498    4.0907    1.7038    3.7166
    1.2896    0.4976    3.9007    2.1549    4.5737
    1.4757    0.3740    3.2786    2.4925    5.4779
    1.4265    0.2978    2.9607    2.3440    3.9290 ];
metric_par.csf_lums = [0.02 0.2 2 20 150];

metric_par.csf_sa = [29.4495 00.188 0.5961 1.6212]; %p5-p8
metric_par.csf_sr_par = [1.1732 1.1478 1.2167 0.5547 2.9899 1.1]; % rod sensitivity function


par = [0.061466549455263 0.99727370023777070]; % old parametrization of MTF
metric_par.mtf_params_a = [par(2)*0.426 par(2)*0.574 (1-par(2))*par(1) (1-par(2))*(1-par(1))];
metric_par.mtf_params_b = [0.028 0.37 37 360];

metric_par.quality_band_freq = [15 7.5 3.75 1.875 0.9375 0.4688 0.2344];

% Using log( X + epsilon );
metric_par.quality_band_w = [0.8299    0.1555    0.0235    0.0488   -0.0992    0.1508   -0.1093];
%metric_par.quality_band_w = [0.7848    0.1231    0.0492    0.0015   -0.0161    0.0849   -0.0274];

%metric_par.quality_band_w = [0.1778    0.2325   -0.0203    0.0975   -0.0170    0.0260    0.0414];
%metric_par.quality_band_w = [137.3990 33.4064 1.5838 5.0754 -11.0824 20.5270 6.0479 ];
metric_par.quality_logistic_q1 = 4.745;
metric_par.quality_logistic_q2 = 0.1755;

metric_par.calibration_date = '18 Jan 2011';

metric_par.surround_l = -1; %use mean image luminance

% process options
i = 1;
while( i <= length( options ) )
    if( strcmp( options{i}, 'pixels_per_degree' ) )
        i = i+1;
        metric_par.pix_per_deg = options{i};
    elseif( strcmp( options{i}, 'viewing_distance' ) )
        i = i+1;
        metric_par.view_dist = options{i};
    elseif( strcmp( options{i}, 'peak_sensitivity' ) )
        i = i+1;
        metric_par.sensitivity_correction = daly_peak_contrast_sens / 10.^(-options{i});
    else
        % all other options
        metric_par.(options{i}) = options{i+1};
        i = i+1;
    end
    i = i+1;
end

end
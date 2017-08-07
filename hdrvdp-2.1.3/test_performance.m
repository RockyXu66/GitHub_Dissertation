% test metric performance

test = pfs_read_luminance( 'threshold_stimuli/hex_patterns/stimuli/t_12.exr' );
reference = pfs_read_luminance( 'threshold_stimuli/hex_patterns/stimuli/r_12.exr' );

tic

options = {}; %{ 'no_masking', 'true' };

res = hdrvdp( test, reference, 'luminance', 30, options );

toc

display( sprintf( 'P_det = %g\nQ_MOS = %g', res.P_det, res.Q_MOS) );

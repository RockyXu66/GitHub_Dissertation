function map = hdrvdp_visualize( P, img, target, colormap )
% HDRVDP_VISUALIZE produces color or grayscale visualization of the metric
% predictions
%
% map = HDRVDP_VISUALIZE( P )
% map = HDRVDP_VISUALIZE( P, img )
% map = HDRVDP_VISUALIZE( P, img, target )
% map = HDRVDP_VISUALIZE( P, img, target, colormap )
%
% The function visualizes error map P and optionally adds a context image
% img. Visual differences P must be within range 0-1. If you do not want to
% show a context image, pass an empty matrix ([]) instead.
%
% 'target' parameter specifies desired output device and can be one of:
%
%   'screen' - (default) - the map be shown on a color screen. The map will
%   contain good reproduction of the context image 'img'.
%  
%   'print' - the map can be printed on a gray-scale printer, so color
%   information will be lost. If this target is selected, a color map used
%   to represent errors 'P' will contain luma differences in addition to
%   color. To ensure that the context image does not interfere with errors,
%   only low-contrast and high frequency content of the image will be
%   preserved. 
%
% 'colormap' parameter can be one of:
%
%   'trichromatic' - this is default if no colormap is specified. Errors are
%   represented as multiple colors: blue, cyan, green, yellow and red, which
%   correspond to P equal 0, 0.25, 0.5, 0.75 and 1.
%
%   'dichromatic' - more appropriate if observes may be color deficient. The
%   hue changes from cyan (0.0) to gray (0.5) and then yellow (1.0). The look
%   of images for color-deficient observers can be tested at:
%   http://www.colblindor.com/coblis-color-blindness-simulator/
%
%   'monochromatic' - use only grayscale. Makes sense only with
%   target='print' or when no context image is specified
% 
%
% Tone-mapping is applied to img to reduce the dynamic range so that highly
% saturated colors can be used also in bright image regions. Tone-mapping
% will enhance contrast if it is lower than the available contrast. This
% behavior is useful for images that have contrast that is near detection
% threshold (e.g. ModelFest stimuli).
%
% The function returns a gamma-corrected sRGB image.
%
% Legend for the color-scales can be found in the color_scales
% directory.
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


if( ~exist( 'colormap', 'var' ) )
    colormap = 'trichromatic';
end

if( ~exist( 'target', 'var' ) )
    target = 'screen';
end

if( ~exist( 'img', 'var' ) || isempty(img) )
    tmo_img = ones(size(P))*0.5;
elseif( strcmp( target, 'print' ) )   
    l = log_luminance( img );
    hp_img = (l - blur_gaussian( l, 2 ) + mean(l(:)));
    tmo_img = vis_tonemap(hp_img, 0.1) + 0.5;
elseif( strcmp( target, 'screen' ) )
    tmo_img = vis_tonemap( log_luminance(img), 0.6 );
else
    error( 'Unknown target: %s', target );
end

    
    if( strcmp( colormap, 'trichromatic' ) || strcmp( colormap, 'print' ) )
        
        color_map = [0.2  0.2  1.0;
            0.2  1.0  1.0;
            0.2  1.0  0.2;
            1.0  1.0  0.2;
            1.0  0.2  0.2];
        
        color_map_in = [0 0.25 0.5 0.75 1];
        
    elseif( strcmp( colormap, 'dichromatic' ) )
        
        color_map = [0.2  1.0  1.0;
            1.0  1.0  1.0
            1.0  1.0  0.2];
        
        color_map_in = [0 0.5 1];

    elseif( strcmp( colormap, 'monochromatic' ) )
        
        color_map = [1.0  1.0  1.0;
            1.0  1.0  1.0];
        
        color_map_in = [0 1];        
        
    else
        error( 'Unknown colormap: %s', colormap );
    end
    

    if( strcmp( target, 'screen' ) )
        color_map_l = color_map * [0.2126 0.7152 0.0722]'; %sum(color_map,2);
        color_map_ch = color_map ./ repmat( color_map_l, [1 3] );
    else
        if( strcmp( colormap, 'monochromatic' ) )
            color_map_l = (color_map * [0.2126 0.7152 0.0722]') ./ color_map_in';
        else
            % luminance map start at 0.3, so that colors are visible
            color_map_l = (color_map * [0.2126 0.7152 0.0722]') ./ (color_map_in'*0.8+0.2);
        end
        color_map_ch = color_map ./ repmat( color_map_l, [1 3] );
    end
    
    %The line belows display pixel values
    %round(min( color_map_ch*255, 255 ))
    
    map = zeros( size(P,1), size(P,2), 3 );
    map(:,:,1) = interp1( color_map_in, color_map_ch(:,1), P );
    map(:,:,2) = interp1( color_map_in, color_map_ch(:,2), P );
    map(:,:,3) = interp1( color_map_in, color_map_ch(:,3), P );
    %map(:,:,3) = 1 - map(:,:,2) - map(:,:,1);
    
    %map = repmat( tmo_img, [1 1 3] );
    map = map .* repmat( tmo_img, [1 1 3] );

end

function l = log_luminance( X )

if( size(X,3) == 3 )
    Y = X(:,:,1) * 0.212656 + X(:,:,2) * 0.715158 + X(:,:,3) * 0.072186;
else
    Y = X;
end

Y(Y<=0) = min(Y(Y>0));
l = log(Y);

end

function tmo_img = vis_tonemap( b, dr )
   
    t = 3;
    
    b_min = min(b(:));
    b_max = max(b(:));
    
    b_scale = linspace( b_min, b_max, 1024 );
    b_p = hist( b(:), b_scale );
    b_p = b_p / sum( b_p(:) );
    
    sum_b_p = sum( b_p.^(1/t) );
    dy = b_p.^(1/t) / sum_b_p;
    
    v = cumsum( dy )*dr + (1-dr)/2;
    
    tmo_img = interp1( b_scale, v, b );
end

function Y = blur_gaussian( X, sigma )
ksize2 = round(sigma*3);
gauss_1d = exp( -(-ksize2:ksize2).^2/(2*sigma^2) );
gauss_1d = gauss_1d/sum(gauss_1d);

Y = conv2( X, gauss_1d, 'same' );
Y = conv2( Y, gauss_1d', 'same' );

end
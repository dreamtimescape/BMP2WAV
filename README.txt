BMP2WAV
Version 1.0 - 2015-09-10
Encode an image (BMP) to an audio file (WAV)
using the spectrum (view with a spectrogram)
-----------------------------------------------------------------
Usage: BMP2WAV -i input.bmp -o output.wav
  -i [File]   Input image (BMP)
  -o [File]   Output audio (WAV 44.1khz 16bit mono)
  -t [float]  Total time duration of output file (s). If omitted, the duration
              will be as long as necessary to make each pixel contain 1.5 cycles
              of the Output low frequency (or use -c to specify value)
  -c [float]  Number of Output Low frequency cycles per pixel - as an alternative
              to specifying the duration, use this option to specify the number of
              cycles of the lowest frequency should output per pixel.
  -p [int]    Samples per Pixel - as an alternative to specifying the duration,
              or number of low frequency cycles, use this option to specify the 
              number of audio samples per pixel.
  -a [int]    Max amplitude per sample (default: 30000) Max: 32768
  -l [int]    Output low frequency (default: 440)
  -h [int]    Output high frequency (default: 19800)
  -x          No silence (encode 0,0,0 with lowest amplitude - not silence)
  -b [string] Background/Transparent color (set to silence) - hex string of RGB
              (ie: 00FF00 => R=0,G=255,B=0)  Note this value will
              replace 0,0,0 as the 'zero amplitude' value, moving 0,0,0
              up to a non-zero brightness.  If -b is not specified, there
              is no background color distinctions made.  Set to 'tl' to use
              the top-left pixel value for the background RGB value.
  -n          Negative (Invert image color)
  -d [float]  dB value for amplitude scaling (default 50). This specifies
              the decibel range between black and white on the color scale.
  -m [int]    Interpolation mode for amplitude between samples
                0 = Hermite (4-point, 3rd order ) (best quality) (default)
                1 = Cubic (4-point, 3rd order)
                2 = Cosine (2-point)
                3 = Linear (2 point, first order)
                4 = None (worst quality)
                5 = Reserved (Test)
                else = Hermite
  -z          Interpolate constant amplitude values - if you provide this 
              option, it will use the interpolation formula even if all the 
              samples are the same amplitude.
  -?          Display this help message
-------------------------------------------------------------------

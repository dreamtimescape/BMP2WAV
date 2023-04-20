#!/usr/bin/perl
#------------------------------------------
#
#Convert bitmap image into audio WAV file that
#can be seen in a spectrogram
#
#------------------------------------------
use Audio::Wav;
use GD;
use Getopt::Std;
use Image::BMP;
use List::Util qw( max min sum );
use Math::Trig;
use POSIX;
use Time::HiRes qw( time );
use strict;
use warnings;

# Elapsed time counter
my $start_time = time();

# Set to nonzero for debug print statements.
my $debug = 0;

#Sample rate for CD quality
#(The only frequency supported as of now)
my $samplerate = 44100;
my $bits_sample = 16;

#Parse command line arguments
my ($input,$output,$duration,$highamp,$lowfreq,$highfreq,$inverse,$bkgnd,$nosilence,$mode,$interpconstamp,$pixelsarg,$lowfcycles,$dbmax) = parseArgs();

print "\n---------------------------\n";
print   " BMP2WAV - Begin\n";
print "---------------------------\n";
print "Input: $input \nOutput: $output \n";
print "Max amplitude per sample: $highamp \n";
print "Image color inverted\n" if ($inverse);
print "dB color scale range: $dbmax\n";

my $modestring = "Hermite";
if ( $mode == 1 ) {
    $modestring = "Cubic";
}
elsif ( $mode == 2 ) {
    $modestring = "Cosine";
}
elsif ( $mode == 3 ) {
    $modestring = "Linear";
}
elsif ( $mode == 4 ) {
    $modestring = "None";
}
elsif ( $mode == 5 ) {
    $modestring = "Test Mode - no calculations";
}
elsif ( $mode == 6 ) {
    $modestring = "4-point, 4th order optimal, 2x";
}
elsif ( $mode == 7 ) {
    $modestring = "4-point, 4th order optimal, 4x";
}
elsif ( $mode == 8 ) {
    $modestring = "4-point, 4th order optimal, 8x";
}
elsif ( $mode == 9 ) {
    $modestring = "4-point, 4th order optimal, 16x";
}
elsif ( $mode == 10 ) {
    $modestring = "4-point, 4th order optimal, 32x";
}
elsif ( $mode == 11 ) {
    $modestring = "Sine Window pixel scaling, zero at edges";
}
elsif ( $mode == 12 ) {
    $modestring = "Hann Window pixel scaling, zero at edges";
}
elsif ( $mode == 13 ) {
    $modestring = "Nutall Window pixel scaling, zero at edges";
}
elsif ( $mode == 14 ) {
    $modestring = "Hann Window half pixel transitions";
}
elsif ( $mode == 15 ) {
    $modestring = "Nutall Window half pixel transitions";
}
else {
    $mode = 0;
    $modestring = "Hermite";
}
print "Interpolation: $modestring\n";

print "Interpolate between constant amplitude samples: ";
if ( $interpconstamp ) {
    print "Yes\n";
}
else {
    print "No\n";
}

if ( $nosilence ) {
    print "No Silence mode: RGB(0,0,0) will be encoded at minimum amplitude\n";
}
if ( lc $bkgnd eq "tl" ) {
    print "Upper Left Pixel will determine Background color.\n";
}

#Open BMP
my $img = new Image::BMP(file=> $input);

#Create file to write to
#open(WAVEFILE, ">$output");

#Get width and height
my $width = $img->{Width}; 
my $height = $img->{Height};
print "Image Width: $width Height: $height \n";

#scale frequency per height;
my $freqrange = $highfreq-$lowfreq;
my $interval = $freqrange / ($height-1);
my @frequencies;

# For 9-pixel tall images, assume we're doing text
if ( $height == 9 ) {
    # If they specify $highfreq lower than $lowfreq, assume they want a 
    # Tone set centered on $lowfreq with $highfreq as the frequency range
    if ( $freqrange < 0 ) {
        $freqrange = $highfreq;
        $interval = $highfreq/($height-1);
        print "9-pixel high font: Using $lowfreq as center tone with $highfreq total range.\n";
        $frequencies[4]=$lowfreq;
        for ( my $i=1; $i<=4; $i++ ) {
            $frequencies[4+$i]=$lowfreq+$i*$interval;
            $frequencies[4-$i]=$lowfreq-$i*$interval;
        }
        $lowfreq = $frequencies[0];
        $highfreq = $frequencies[8];
    }
    # If they specify $highfreq exactly equal to $lowfreq, use a chord
    elsif ( $freqrange = 0 ) {
        print "Using A chord for 9-pixel high input\n";
        # A major = A,C#,E
        # A minor = A,C,E
        $frequencies[0] = 440;      # A4
        $frequencies[1] = 659.3;    # E5
        $frequencies[2] = 880;      # A5
        $frequencies[3] = 1047;     # 1047=C6, 1109=C#6
        $frequencies[4] = 1319;     # E6
        $frequencies[5] = 1760;     # A6
        $frequencies[6] = 2093;     # 2093=C7, 2217=C#7
        $frequencies[7] = 2637;     # E7
        $frequencies[8] = 3520;     # A7
    }
    else {
        for ( my $i=0; $i<$height; $i++ ) {
            $frequencies[$i] = $lowfreq + $i*$interval;
        }
    }
}
else {
    for ( my $i=0; $i<$height; $i++ ) {
        $frequencies[$i] = $lowfreq + $i*$interval;
    }
}

print "Low Freq: $lowfreq Hz, High Freq: $highfreq Hz\n";
printf ( "Frequency Interval: %.3f Hz\n", $interval);

# Set the low and high limits for outputting pixels
# Total pixels output = $highlimit - $lowlimit
my $lowlimit = 1;
my $highlimit = $width+1;
# For half-pixel transition window functions, need extra pixels at beginning and end since edge values are fractional
if ( ( $mode == 14 ) || ( $mode == 15 ) ) {
    $lowlimit = 0;
    $highlimit = $width+2;
}
my $outputwidth = $highlimit - $lowlimit;

# Default for samples/pixel.
my $samplespixel = 400;

# Calculate duration if not specified
if ( $duration == 0 ) {
    # If samples per pixel was provided, use that number
    if ( $pixelsarg == 0 ) {
        # If low frequency cycles is provided, use that value, otherwise zero it out.
        if ( $lowfcycles ) {
            print "Low Frequency Cycles Mode\n";
        }
        else {
            print "Default Duration Mode\n";
            $lowfcycles = 1.5;
        }
        print "$lowfcycles cycles per pixel for $lowfreq Hz.\n";
        # Each pixel to contain $lowfcycles cycles of the lowest frequency
        $duration = ($width * $lowfcycles * ($samplerate/$lowfreq))/$samplerate;
        $samplespixel =  POSIX::ceil($samplerate * $duration/ ($outputwidth));  
    }
    else {
        print "Pixelsize mode\n";
        $samplespixel = $pixelsarg;
    }
}
else {
    print "Specified Duration Mode.\n";
    # Calculate the samples per pixel based on image size and total duration.
    $samplespixel =  POSIX::ceil($samplerate * $duration/ ($outputwidth));  
}
# Recalculate duration based on samples per pixel due to ceiling function
$duration = sprintf("%.3f",($samplespixel * ($outputwidth))/$samplerate);

printf ( "Total duration: %.3f sec\n", $duration );
print "Samples per Pixel: $samplespixel \n";


# Sine Look up table
my @sinelut;
my $lutsize = 2048;
#print "Calculating Sine Look up table\n";
for ( my $i=0;$i<$lutsize;$i++ ) {
    $sinelut[$i] = sin( $i*(2*pi)/($lutsize) );
}
#print "Done\n";

my $checkbright = 0;

my $maxbright=0;
my $minbright=$highamp;
#Get max/min brightness
if ( $checkbright ) {
    my @brightness;
    for (my $x=0;$x<$width;$x++) {
        for (my $y=0;$y<$height;$y++) {
            $brightness[$x][$y]= color($img->xy_rgb($x,$y));
            if ( $debug ) {
                print "Brightness (X,Y): ($x,$y): $brightness[$x][$y]\n";
            }
            if ( $brightness[$x][$y] > $maxbright ) {
                $maxbright = $brightness[$x][$y];
            }
            if ( $brightness[$x][$y] < $minbright ) {
                $minbright = $brightness[$x][$y];
            }
        }
    }
    print "Image Brightness Min: $minbright Max: $maxbright\n";
    if ( $debug ) {
        printf ("Brightness (0,0,0) = %f\n", color( 0, 0, 0) );
        printf ("Brightness (1,0,0) = %f\n", color( 1, 0, 0) );
        printf ("Brightness (0,1,0) = %f\n", color( 0, 1, 0) );
        printf ("Brightness (0,0,1) = %f\n", color( 0, 0, 1) );
        printf ("Brightness (1,1,1) = %f\n", color( 1, 1, 1) );
    }
}


#This is the main algorithm
my @audioOut;
my $cursample= 0;
my @phaseshift;
my @prevamp;
my @zerowave = (0) x $samplespixel;
my @phasedeltas;
my @maxamps;
my @freqcount;
my $maxfreqs = 0;

# Set pure to nonzero for actual sine calculations.  If zero, look-up table is used.
my $pure = 1;

# Null value (minimum amplitude) for pre and post image pixels
# Black for normal images
my $null = color(0,0,0);
if ( $inverse ) {
    # White for inverted
    $null = color(255,255,255);
}

#Move across the file one pixel at a time
# Every pixel row is a different oscillator
# Every pixel column is a sample set of weighted sums of those oscillators
for (my $x=0;$x<$width+2;$x++) {
    # Keep track of current sample for phase
    $cursample = ($x) * $samplespixel;

    #Calculate rows;
    my @audiorow;
    my $pos=0;
    $maxamps[$x] = 0;
    for (my $y=0;$y<$height;$y++) {
        
        # Calculate this row's random phase shift
        if ( $pure ) {
            $phaseshift[$y] = rand(2*pi) unless defined $phaseshift[$y];
        }
        else {
            $phaseshift[$y] = rand($lutsize-1) unless defined $phaseshift[$y];
        }

        my $amp=0;
        my $prev=0;
        # Check for top-left background color condition
        if (($x+$y==0) && ( lc $bkgnd eq "tl" )) {
            my ($r,$g,$b) = $img->xy_rgb($x,$y);
            $bkgnd = sprintf("%02x",$r) . sprintf("%02x",$g) . sprintf("%02x",$b);
            $bkgnd = uc $bkgnd;
            print "Setting Background to: $bkgnd ($r,$g,$b)\n";
            # Redo the null value
            $null = color( $r, $g, $b );
        }
        my $yinv = $height - $y - 1;
        #Generate a sine wave 
        #my $freq=$yinv*$interval + $lowfreq;
        my $freq=$frequencies[$yinv];

        # Incremental phase for this frequency.
        $phasedeltas[$y] = ($freq*2*pi)/$samplerate;

        if ( $x < $width ) {
            #Get the amplitude of this pixel
            $amp = color($img->xy_rgb($x,$y));
        }
        else {
            # Set the post-image pixels to match the last pixel.
            #$amp = $prevamp[$y][3];
            # Set the post image pixels to null value
            $amp = $null;
        }
        if ( $x == 0 ) {
            # Set the pre-image pixels to match the first pixel.
            #@{ $prevamp[$y] } = ( $amp ) x 4;
            # Set the pre-image pixels to null value
            @{ $prevamp[$y] } = ( $null ) x 4;
        }

        # Shift the previous value array
        shift @{ $prevamp[$y] };
        # Add the new amplitude 
        push @{ $prevamp[$y] }, $amp;
        
        # Read in 2 values before calculating first pixel's output
        # We'll end up with an extra transition pixel before the first
        if ( $x >= $lowlimit ) {
            # Skip if all amplitude values are zero.
            if ( sum( @{ $prevamp[$y] } ) ) {
                $audiorow[$pos++] = genSine( $phasedeltas[$y], $phaseshift[$y], $cursample, $samplespixel, $samplerate, @{$prevamp[$y]} );
            }
        }
        $maxamps[$x] += $amp;
        $freqcount[$x] += ( $amp != 0 );
        
    }
    # Start outputting data on $lowlimit column [0-based] and output until $highlimit (handles 'null' pre- and post- image pixels and window functions)
    if ( ( $x >= $lowlimit ) && ( $x <= $highlimit ) ) {
    
        # If there was nothing to output, create a 'zero' wave
        if ($pos==0) {
            $audiorow[$pos++] = \@zerowave;
            #$audiorow[$pos++] = genSine( $lowfreq, 0, 0, $samplespixel, $samplerate, ((0) x 4) );
        }

        #Add all the samples for the column
        for(my $i=0;$i<$samplespixel;$i++) {
            my $index = $i + ($x-$lowlimit) * $samplespixel;
            for(my $j=0;$j<@audiorow;$j++) {
                $audioOut[$index] += $audiorow[$j][$i];
                #print "index: $index, j:$j, i:$i\n";
            }
        }
        
        # Keep track of how many frequencies were summed in the most-summed column
        if ( $freqcount[$x] > $maxfreqs ) {
            $maxfreqs = $freqcount[$x];
        }
    }

}

# Scale the final output - slightly less than peak sample value as it's not likely the very top of the sine wave.
my $maxval = max( @audioOut ) * 1.05;
print "Max actual amplitude: $maxval\n";

print "Max summed freqs: $maxfreqs\n";
# Worst case amplitude - if all sine waves in a column were summed with constructive interference
# this would be the amplitude.  Need to scale it down slightly due to random phase difference?
# What value? What's the amplitude of a sum of randomly phased known amplitude sinusoids?  That's a complex question.
my $maxworst = max( @maxamps );
print "Max theoretical amp: $maxworst\n" ;
# Let's just say if the maximum theoretical is more than T times the actual maximum, 
# just use the actual maximum * T.
# Pretty much making up values for T.
my $Tval = sqrt(2);
if ( $maxworst > $Tval * $maxval ) {
    print "Theoretical max more than $Tval times actual, using actual*$Tval\n";
    $maxworst = $Tval * $maxval;
}
# If something went wrong ( I think I saw it happen once) make sure the actual highest known amplitude is used.
elsif ( $maxval > $maxworst ) {
    print "WHOAH!! Maximum Value Exceeds Theoretical Maximum - Check your theory!\n";
    $maxworst = $maxval;
}
print "Scaling $maxworst to Maximum amplitude of $highamp\n";

if ( $maxworst ) {
    #foreach my $x (@audioOut) { $x *= $highamp/$maxval };
    foreach my $x (@audioOut) { $x *= $highamp/$maxworst };
}


print "\nGenerating wave file \n";
#Write samples to wave file
#print WAVEFILE SimpleWave::genWave(\@audioOut);
my $wav = Audio::Wav->new;
my $write = $wav->write( $output,
                        {
                        bits_sample => $bits_sample,
                        sample_rate => $samplerate,
                        channels    => 1,
                        }
                       );
#Write samples to wave file
# 256K at a time
my $chunksize = 262144;
my $chunks = int(@audioOut/$chunksize);
my $remains = @audioOut - $chunks*$chunksize;
print "Splitting into $chunks chunks of $chunksize plus $remains\n";
for ( my $i=0;$i<$chunks;$i++ ) {
    $write->write( @audioOut[($i*$chunksize) .. ($i*$chunksize + ($chunksize - 1))] );
}
if ( $remains ) {
    $write->write( @audioOut[($chunks*$chunksize) .. ($chunks*$chunksize + ($remains - 1))] );
}
$write->finish;

my $end_time = time();
print "\n---------------------------\n";
print " BMP2WAV - Complete\n";
printf ("              Time: %.3f\n", $end_time - $start_time );
print "---------------------------\n";

#-----------------------------------------
#
#Generate a sine wave for the given 
#frequency and amplitude 
#
# Previous amplitude provided to keep
# sine wave at that amplitude until the
# first zero-crossing, then switch to the
# new amplitude.
#
#-----------------------------------------
sub genSine {

    #my ($frequency, $phase, $cursample, $samples, $samplerate, @amplitudes) = @_;
    my ($phasedelta, $phase, $cursample, $samples, $samplerate, @amplitudes) = @_;
    
    my @audioArray;
    my $prevamp = 0;
    my $amplitude = 0;
    my $nextamp = 0;
    my $slope = 0;
#    my $foverfc = 0;
    my $totphase = 0;
#    my $phasedelta = 0;
    my $y = 0;
    my $transamp=0;
    my $ampphase = 0;
    my $ampphasedelta=0;
    my $mymode = $mode;

    # Skip the whole shebang if all amplitudes are zero
    if ( sum( @amplitudes ) ) {
        $prevamp = $amplitudes[1];
        $amplitude = $amplitudes[2];
        $nextamp = $amplitudes[3];
        $slope = ($amplitude - $prevamp)/$samples;
#        $foverfc = ($frequency/$samplerate);
        $totphase = 0;
#        $phasedelta = $foverfc*2*pi;
        $y = 0;
        $transamp=$prevamp;
        $ampphase = 0;
        $ampphasedelta=0;
        
        if ( $pure ) {
            # If the previous and next amplitudes are the same, Don't interpolate unless interpconstamp is set
            if ( (!$interpconstamp) && ($prevamp == $amplitude ) ) {
                # For the window functions, need to interpolate every sample
                if ( $mode <= 10 ) {
                    # No interpolation.
                    $mymode = 4;
                }
            }
            
            # MODE: 0 = Hermite
            #       1 = Cubic
            #       2 = Cosine
            #       3 = Linear
            #       4 = None
            #       5 = Testing - no calculations
            #       6 = Hermite with half pixel transitions
            #    else = Hermite
            if ( $mymode == 0 ) {
                # Also covers $mode == 0
                # Hermite interpolation 4-point, 3rd order
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $c0 = $amplitudes[1];
                my $c1 = ($amplitudes[2] - $amplitudes[0])/2;
                my $c2 = $amplitudes[0] - (5 * $amplitudes[1]/2) + ( 2 * $amplitudes[2] ) - ( $amplitudes[3]/2 );
                my $c3 = ( ($amplitudes[3] - $amplitudes[0])/2 ) + ( 3 * ( $amplitudes[1] - $amplitudes[2] )/2 );
                my $tdelt = 1/$samples;
                my $t = 0;
                for(my $i=0;$i<$samples;$i++) {
                    # Interpolated amplitude
                    $transamp = ((((($c3*$t) + $c2) * $t) + $c1) * $t) + $c0;
                    # Post Increment interval
                    $t += $tdelt;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            elsif ( $mymode == 1 ) {
                # Cubic interpolation
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $a0 = $amplitudes[3] - $amplitudes[2] - $amplitudes[0] + $amplitudes[1];
                my $a1 = $amplitudes[0] - $amplitudes[1] - $a0;
                my $a2 = $amplitudes[2] - $amplitudes[0];
                my $a3 = $amplitudes[1];
                for(my $i=0;$i<$samples;$i++) {
                    # Initialize to zero.
                    $y = 0;
                    # Interpolated amplitude
                    my $t = $i/$samples;
                    $transamp = ($a0 * ($t*$t*$t)) + ($a1 * ($t*$t)) + ($a2 * $t) + $a3;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ( $mymode == 2 ) {
                # Shape amplitude transition with cosine curve.
                $slope = ($amplitude - $prevamp)/2;
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                $transamp=$prevamp + $slope;
                $ampphase = 0;
                $ampphasedelta = ( pi/$samples );
                for(my $i=0;$i<$samples;$i++) {
                    # Initialize to zero.
                    $y = 0;
                    # Scaled amplitude
                    $ampphase += $ampphasedelta;
                    $transamp = $prevamp + $slope*(1-cos($ampphase));
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ( $mymode == 3 ) {
                # Linear interpolation
                $slope = ($amplitude - $prevamp)/$samples;
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                $transamp=$prevamp;
                for(my $i=0;$i<$samples;$i++) {
                    # Scaled amplitude
                    $transamp += $slope;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            elsif ( $mymode == 4 ) {
                # No interpolation
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                $transamp=$amplitude;
                for(my $i=0;$i<$samples;$i++) {
                    # Increment phase
                    $totphase += $phasedelta;
                    # Calculate the sample value.
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            elsif ( $mymode == 5 ) {
                # Test with no calculations - doesn't produce a sine wave.
                for(my $i=0;$i<$samples;$i++) {
                    $y = $amplitude;
                    $audioArray[$i] = $y;
                }
            }
            elsif ($mymode == 6) {
                # Optimal 2x (4-point, 4th order)
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $even1 = $amplitudes[2] + $amplitudes[1];
                my $odd1 =  $amplitudes[2] - $amplitudes[1];
                my $even2 = $amplitudes[3] + $amplitudes[0];
                my $odd2 =  $amplitudes[3] - $amplitudes[0];
                my $c0 = $even1*0.45645918406487612 + $even2*0.04354173901996461;
                my $c1 = $odd1*0.47236675362442071 + $odd2*0.17686613581136501;
                my $c2 = $even1*-0.253674794204558521 + $even2*0.25371918651882464;
                my $c3 = $odd1*-0.37917091811631082 + $odd2*0.11952965967158000;
                my $c4 = $even1*0.04252164479749607 + $even2*-0.04289144034653719;
                for(my $i=0;$i<$samples;$i++) {
                    my $x = $i/$samples;
                    my $z = $x - 0.5;
                    $transamp = ((($c4*$z + $c3) * $z + $c2) * $z + $c1) * $z + $c0;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ($mymode == 7) {
                # Optimal 4x (4-point, 4th order)
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $even1 = $amplitudes[2] + $amplitudes[1];
                my $odd1 =  $amplitudes[2] - $amplitudes[1];
                my $even2 = $amplitudes[3] + $amplitudes[0];
                my $odd2 =  $amplitudes[3] - $amplitudes[0];
                my $c0 = $even1*0.46567255120778489 + $even2*0.03432729708429672;
                my $c1 = $odd1*0.53743830753560162 + $odd2*0.15429462557307461;
                my $c2 = $even1*-0.251942101340217441 + $even2*0.25194744935939062;
                my $c3 = $odd1*-0.46896069955075126 + $odd2*0.15578800670302476;
                my $c4 = $even1*0.00986988334359864 + $even2*-0.00989340017126506;
                for(my $i=0;$i<$samples;$i++) {
                    my $x = $i/$samples;
                    my $z = $x - 0.5;
                    $transamp = ((($c4*$z + $c3) * $z + $c2) * $z + $c1) * $z + $c0;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ($mymode == 8) {
                # Optimal 8x (4-point, 4th order)
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $even1 = $amplitudes[2] + $amplitudes[1];
                my $odd1 =  $amplitudes[2] - $amplitudes[1];
                my $even2 = $amplitudes[3] + $amplitudes[0];
                my $odd2 =  $amplitudes[3] - $amplitudes[0];
                my $c0 = $even1*0.46771532012068961 + $even2*0.03228466824404497;
                my $c1 = $odd1*0.55448654344364423 + $odd2*0.14851181120641987;
                my $c2 = $even1*-0.250587283698110121 + $even2*0.25058765188457821;
                my $c3 = $odd1*-0.49209020939096676 + $odd2*0.16399414834151946;
                my $c4 = $even1*0.00255074537015887 + $even2*-0.00255226912537286;
                for(my $i=0;$i<$samples;$i++) {
                    my $x = $i/$samples;
                    my $z = $x - 0.5;
                    $transamp = ((($c4*$z + $c3) * $z + $c2) * $z + $c1) * $z + $c0;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ($mymode == 9) {
                # Optimal 16x (4-point, 4th order)
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $even1 = $amplitudes[2] + $amplitudes[1];
                my $odd1 =  $amplitudes[2] - $amplitudes[1];
                my $even2 = $amplitudes[3] + $amplitudes[0];
                my $odd2 =  $amplitudes[3] - $amplitudes[0];
                my $c0 = $even1*0.46822774170144532 + $even2*0.03177225758005808;
                my $c1 = $odd1*0.55890365706150436 + $odd2*0.14703258836343669;
                my $c2 = $even1*-0.250153411893796031 + $even2*0.25015343462990891;
                my $c3 = $odd1*-0.49800710906733769 + $odd2*0.16600005174304033;
                my $c4 = $even1*0.00064264050033187 + $even2*-0.00064273459469381;
                for(my $i=0;$i<$samples;$i++) {
                    my $x = $i/$samples;
                    my $z = $x - 0.5;
                    $transamp = ((($c4*$z + $c3) * $z + $c2) * $z + $c1) * $z + $c0;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ($mymode == 10) {
                # Optimal 32x (4-point, 4th order)
#                $foverfc = ($frequency/$samplerate);
#                $phasedelta = $foverfc*2*pi;
                $totphase = $phasedelta*$cursample;
                my $even1 = $amplitudes[2] + $amplitudes[1];
                my $odd1 =  $amplitudes[2] - $amplitudes[1];
                my $even2 = $amplitudes[3] + $amplitudes[0];
                my $odd2 =  $amplitudes[3] - $amplitudes[0];
                my $c0 = $even1*0.46835497211269561 + $even2*0.03164502784253309;
                my $c1 = $odd1*0.56001293337091440 + $odd2*0.14666238593949288;
                my $c2 = $even1*-0.250038759826233691 + $even2*0.25003876124297131;
                my $c3 = $odd1*-0.49949850957839148 + $odd2*0.16649935475113800;
                my $c4 = $even1*0.00016095224137360 + $even2*-0.00016095810460478;
                for(my $i=0;$i<$samples;$i++) {
                    my $x = $i/$samples;
                    my $z = $x - 0.5;
                    $transamp = ((($c4*$z + $c3) * $z + $c2) * $z + $c1) * $z + $c0;
                    # Calculate the sample value.
                    $totphase += $phasedelta;
                    $y = $transamp*sin($totphase + $phase );
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $y;
                }
            }
            elsif ( $mymode == 11 ) {
                # Sine Window
                # Scale new amplitude by sin(x) where x ranges from 0 to pi over the pixel.
                $totphase = $phasedelta*$cursample;
                $ampphasedelta = ( pi/$samples );
                $ampphase = $ampphasedelta * $cursample;
                for(my $i=0;$i<$samples;$i++) {
                    # Increment amplitude phase counter
                    $ampphase += $ampphasedelta;
                    $transamp = $amplitude*sin($ampphase);
                    # Increment phase tracker
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            elsif ( $mymode == 12 ) {
                # Hann Window
                # f = 1/2*(1 - cos(2*pi*n/(N-1)))
                $totphase = $phasedelta*$cursample;
                $ampphase = 0;
                $ampphasedelta = ( 2*pi/($samples-1) );
                my $alpha = 1/2;
                my $beta = 1/2;
                for(my $i=0;$i<$samples;$i++) {
                    # Increment window function phase counter
                    $ampphase += $ampphasedelta;
                    # Calculate Window amplitude * input amplitude
                    $transamp = $amplitude*( $alpha - $beta*(cos($ampphase)) );
                    # Increment phase tracker
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            elsif ( $mymode == 13 ) {
                # Nutall Window
                # theta = (2*pi*n/(N-1));
                # f = a0 - a1*cos(theta) + a2*cos(2*theta) + a3*cos(3*theta);
                # a0= 0.355768; a1= 0.487396; a2= 0.144232; a3= 0.012604;
                $totphase = $phasedelta*$cursample;
                $ampphase = 0;
                $ampphasedelta = ( 2*pi/($samples-1) );
                my $a0= 0.355768;
                my $a1= 0.487396;
                my $a2= 0.144232;
                my $a3= 0.012604;
                for(my $i=0;$i<$samples;$i++) {
                    # Calculate Window amplitude * input amplitude
                    $transamp = $amplitude*( $a0 - $a1*(cos($ampphase)) + $a2*(cos(2*$ampphase)) - $a3*(cos(3*$ampphase)));
                    # Increment phase tracker
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                    # Increment window function phase counter
                    $ampphase += $ampphasedelta;
                }
            }
            elsif ( $mymode == 14 ) {
                # Hann Window half pixel transitions
                # f = 1/2*(1 - cos(2*pi*n/(N-1)))
                $totphase = $phasedelta*$cursample;
                $ampphase = pi/2;
                $ampphasedelta = ( 2*pi/($samples-1) );
                my $alpha = 1/2;
                my $beta = 1/2;
                my $slope1 = ($amplitude - $prevamp);
                my $slope2 = ($nextamp - $amplitude);
                for(my $i=0;$i<$samples;$i++) {
                    if ( $i < $samples/4 ) {
                        # Finish Transition from $prevamp to $amplitude in first quarter
                        $transamp = $prevamp + $slope1*( $alpha - $beta*(cos($ampphase)) );
                        # Increment window function phase counter
                        $ampphase += $ampphasedelta;
                    }
                    elsif ( $i > $samples - $samples/4 ) {
                        # Start Transition from $amplitude to $nextamp in last quarter
                        $transamp = $amplitude + $slope2*( $alpha - $beta*(cos($ampphase)) );
                        # Increment window function phase counter
                        $ampphase += $ampphasedelta;
                    }
                    else {
                        # Reset Window Phase counter for last quarter
                        $ampphase = 0;
                        # Hold at $amplitude for middle half
                        $transamp = $amplitude;
                    }
                    # Increment phase tracker
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            elsif ( $mymode == 15 ) {
                # Nutall Window half pixel transitions
                # theta = (2*pi*n/(N-1));
                # f = a0 - a1*cos(theta) + a2*cos(2*theta) + a3*cos(3*theta);
                # a0= 0.355768; a1= 0.487396; a2= 0.144232; a3= 0.012604;
                $totphase = $phasedelta*$cursample;
                $ampphase = pi/2;
                $ampphasedelta = ( 2*pi/($samples-1) );
                my $a0= 0.355768;
                my $a1= 0.487396;
                my $a2= 0.144232;
                my $a3= 0.012604;
                my $slope1 = ($amplitude - $prevamp);
                my $slope2 = ($nextamp - $amplitude);
                for(my $i=0;$i<$samples;$i++) {
                    if ( $i < $samples/4 ) {
                        # Finish Transition from $prevamp to $amplitude in first quarter
                        $transamp = $prevamp + $slope1*( $a0 - $a1*(cos($ampphase)) + $a2*(cos(2*$ampphase)) - $a3*(cos(3*$ampphase)));
                        # Increment window function phase counter
                        $ampphase += $ampphasedelta;
                    }
                    elsif ( $i > $samples - $samples/4 ) {
                        # Start Transition from $amplitude to $nextamp in last quarter
                        $transamp = $amplitude + $slope2*( $a0 - $a1*(cos($ampphase)) + $a2*(cos(2*$ampphase)) - $a3*(cos(3*$ampphase)));
                        # Increment window function phase counter
                        $ampphase += $ampphasedelta;
                    }
                    else {
                        # Reset Window Phase counter for last quarter
                        $ampphase = 0;
                        # Hold at $amplitude for middle half
                        $transamp = $amplitude;
                    }
                    # Increment phase tracker
                    $totphase += $phasedelta;
                    # Add the scaled amplitude sample to the array.
                    $audioArray[$i] = $transamp*sin($totphase + $phase );
                }
            }
            else {
                # Invalid - main code should prevent this
                print "Invalid Mode!\n";
            }
        }
        else {
            # Attempt at look-up table approach to speed things up.
            my $lutindex = 0;
            my $lutindex_int = 0;
            my $lutindex_dec = 0;
            # Hermite interpolation 4-point, 3rd order
#            $foverfc = ($frequency/$samplerate);
#            $phasedelta = $foverfc*2*pi;
            my $foverfc = $phasedelta/(2*pi);
            $totphase = $phasedelta*$cursample;
            my $c0 = $amplitudes[1];
            my $c1 = ($amplitudes[2] - $amplitudes[0])/2;
            my $c2 = $amplitudes[0] - (2.5 * $amplitudes[1]) + ( 2 * $amplitudes[2] ) - ( $amplitudes[3]/2 );
            my $c3 = ( ($amplitudes[3] - $amplitudes[0])/2 ) + ( 1.5 * ( $amplitudes[1] - $amplitudes[2] ) );
            for (my $i=0;$i<$samples;$i++) {
                $lutindex = (($foverfc)*($i+$cursample)+$phase);
                #print "TOTP: $totphase\n";
                $lutindex = ($lutindex - int($lutindex)) * ($lutsize);
                #print "LUTIND: $lutindex\n";
                $lutindex_int = int($lutindex);
                $lutindex_dec = $lutindex - $lutindex_int;
                if ( 0 ) {
                    $y = $sinelut[$lutindex_int];
                }
                else {
                    $y = ($sinelut[($lutindex_int)] + $lutindex_dec*($sinelut[($lutindex_int+1)&($lutsize-1)]-$sinelut[$lutindex_int]));
                }
                # Interpolated amplitude
                my $t = $i/$samples;
                $transamp = ((((($c3*$t) + $c2) * $t) + $c1) * $t) + $c0;
                $audioArray[$i] = $transamp*$y;
            }
        }
    }
    else {
        # No non-zero amplitude values
        @audioArray = (0) x $samples;
    }

    return \@audioArray;
}

#-----------------------------------------
#
#Convert the RGB value to a amplitude
#
#-----------------------------------------
sub color {
    my ($r,$g,$b) = @_;
    my $total = 0;
    my $amp = $total;
    my $mymax = 100;
    my $CIEL = 1;
    
    # Create string for this color
    my $rgbstr = sprintf("%02x",$r) . sprintf("%02x",$g) . sprintf("%02x",$b);
    
    if ( lc $rgbstr eq lc $bkgnd ) {
        # If the color is the background color, make this the zero value
        # 0,0,0 for regular, 255,255,255 for inverted.
        if ( $inverse ) {
            $r = 255;
            $g = 255;
            $b = 255;
        }
        else {
            $r = 0;
            $g = 0;
            $b = 0;
        }
    }
    elsif ( ($inverse == 0) && (lc $rgbstr eq "000000") ) {
        # Normal scale, the color is black, and it's not the background color
        # Bump it up to non-zero if nosilence is set
        if ( ($nosilence == 1) || ( length( $bkgnd ) > 0 ) ) {
            $r = 0;
            $g = 0;
            $b = 1;
        }
    }
    elsif ( ($inverse == 1) && (lc $rgbstr eq "ffffff") ) {
        # Inverse scale, the color is white, and it's not the background color
        # Bump it down to non-zero nosilence is set
        if ( ($nosilence == 1) || ( length( $bkgnd ) > 0 ) ) {
            $r = 255;
            $g = 255;
            $b = 254;
        }
    }
    if ( $CIEL ) {
        # 100 is max for CIE L* brightness calculation
        $mymax = 100;
        my $var_R = ( $r / 255 );
        my $var_G = ( $g / 255 );
        my $var_B = ( $b / 255 );

        if ( $var_R > 0.04045 ) {
            $var_R = ( ( $var_R + 0.055 ) / 1.055 ) ** 2.4;
        }
        else {
            $var_R = $var_R / 12.92;
        }
        if ( $var_G > 0.04045 ) {
            $var_G = ( ( $var_G + 0.055 ) / 1.055 ) ** 2.4;
        }
        else {
            $var_G = $var_G / 12.92;
        }
        if ( $var_B > 0.04045 ) {
            $var_B = ( ( $var_B + 0.055 ) / 1.055 ) ** 2.4;
        }
        else {
            $var_B = $var_B / 12.92;
        }
        $var_R = $var_R * 100;
        $var_G = $var_G * 100;
        $var_B = $var_B * 100;

        my $X = $var_R * 0.4124 + $var_G * 0.3576 + $var_B * 0.1805;
        my $Y = $var_R * 0.2126 + $var_G * 0.7152 + $var_B * 0.0722;
        my $Z = $var_R * 0.0193 + $var_G * 0.1192 + $var_B * 0.9505;

        my $ref_X = 95.047;
        my $ref_Y = 100.000;
        my $ref_Z = 108.883;

        my $var_X = $X / $ref_X;
        my $var_Y = $Y / $ref_Y;
        my $var_Z = $Z / $ref_Z;
        
        my $CIE_Lstar;
        
        if ( $var_Y == 0 ) {
            $CIE_Lstar = 0;
        }
        else {
            if ( $var_Y > 0.008856 ) {
                $var_Y = $var_Y ** (1/3);
            }
            else {
                $var_Y = ( 7.787 * $var_Y ) + ( 16/116 );
            }

            $CIE_Lstar = ( 116 * $var_Y ) - 16;
        }

        $total = $CIE_Lstar;
    }
    else {
        $mymax = 255;
        $total = 0.2126*$r + 0.7152*$g + 0.0722*$b;
        #$total = ($r + $g + $b)/3;
        #$total = max( $r, $g, $b );
        #$total = max(  $r, $g, $b,  )/2 + min( $r, $g, $b, )/2;
        #$total = 0.30*$r + 0.59*$g + 0.11*$b;
        #print "R: $r, G: $g, B: $b MAX: $total\n";
    }

    #Invert the color if -I is set
    $total = $mymax - $total if($inverse);
    #$amp = 20*log(($total+1))/log(10);
    
    if ( $total == 0 ) {
        $amp = $total;
    }
    else {
        # Get linear amplitude from dB brightness value.
        $amp = (10**(($dbmax*($total))/(20*($mymax))));
#        $amp -= 1;
    }

    return $amp;
}


#-----------------------------------------
#
#Read arguments from the command line;
#
#-----------------------------------------
sub parseArgs {

    #Default arguments
    my $duration = 0;
    my $highamp = 32768;
    my $inverse = 0;
    my $lowfreq = 440;
    my $highfreq = 19800;
    my $bkgnd = "";
    my $nosilence = 0;
    my $mode = 0;
    my $interpconstamp = 0;
    my $pixelsarg = 0;
    my $lowfcycles = 0;
    my $dbmax = 50;

    #Set options for reading
    my %args = ();
    getopts("i:o:t:a:l:h:b:m:p:c:d:nxz\?",\%args);

    help() if($args{"\?"}); 
    help("Need input filename") if(!$args{i}); 
    help("Need output filename") if(!$args{o}); 
    $duration = $args{t} if($args{t}); 
    $highamp = $args{a} if($args{a}); 
    $lowfreq = $args{l} if ($args{l});
    $highfreq = $args{h} if ($args{h});
    $bkgnd = $args{b} if ($args{b});
    $inverse = 1 if($args{n}); 
    $nosilence = 1 if($args{x}); 
    $mode = $args{m} if ($args{m});
    $interpconstamp = 1 if ($args{z});
    $pixelsarg = $args{p} if ($args{p});
    $lowfcycles = $args{c} if ($args{c});
    $dbmax = $args{d} if ($args{d});

    return ($args{i},$args{o},$duration,$highamp,$lowfreq,$highfreq,$inverse,$bkgnd,$nosilence,$mode,$interpconstamp,$pixelsarg,$lowfcycles,$dbmax);
}


#-----------------------------------------
#
#Display help message for program
#
#-----------------------------------------
sub help {

    my ($message) = @_;
    print <<EOF;
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
EOF
    print "Error: $message \n" if($message ne '');
    exit(1);
}


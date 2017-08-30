# Can be called by 
#
#   for i in {1..5}; do index=`printf "%04i" "$i"`; echo "Converting $index..."; python convert_dump_to_input.py v18/run$index/ip1.txt v20/11sec_dist/$i.txt 1000 451209.019545066752471; done;

import sys

DUMP_file = sys.argv[1]
output_file = sys.argv[2]
nbr_particles = int(sys.argv[3])
energy = float(sys.argv[4]) # in MeV

xcorr = -1.999999996185254147462700000000000; # all units 1e-3
pxcorr = 0.000000000150460499225583966566850;
ycorr = -0.000000006411949440904753153327300;
pycorr = -0.169999999800309692377100000000000;

with open(DUMP_file, 'r') as i_f:
    with open(output_file, 'w') as o_f:
        # Reading the last lines -- this might work badly for very large DUMP files
        count = 0
        for line in reversed(i_f.readlines()):
            cols = map(float, line.split())
            x = (cols[3] - xcorr)*1e-3
            px = (cols[4] - pxcorr)*1e-3
            y = (cols[5] - ycorr)*1e-3
            py = (cols[6] - pycorr)*1e-3
            z = cols[7]
            e = (cols[8] + 1.0)*energy
            o_f.write("{0} {1} {2} {3} {4} {5}\n".format(x, px, y, py, z, e))
            #o_f.write("{0} {1} {2} {3} {4} {5}\n".format(0, 0, 0, 0, z, e))
            
            count += 1
            if count == nbr_particles: break

test:
	gunzip -c methyloxirane.fchk.gz > methyloxirane.fchk
	../../polar.x -fchk methyloxirane.fchk -wmin 0.0 -wmax 4000 -NP 10000 -wg 10.0 -vdat test.dat 
	zgrep 'Rot. str.' methyloxirane.out.gz
	zgrep 'IR Inten'  methyloxirane.out.gz
	rm -rf methyloxirane.fchk

icm:
	gunzip -c methyloxirane.fchk.gz > methyloxirane.fchk
	../../polar.x -fchk methyloxirane.fchk -icm_IR   methyloxirane-icm_IR.dat
	../../polar.x -fchk methyloxirane.fchk -icm_VCD  methyloxirane-icm_VCD.dat
	../../polar.x -fchk methyloxirane.fchk -icm_MAG  methyloxirane-icm_MAG.dat
	rm -rf methyloxirane.fchk

robust:
	gunzip -c methyloxirane.fchk.gz > methyloxirane.fchk
	../../polar.x -r -fchk methyloxirane.fchk 
	rm -rf methyloxirane.fchk

all: test icm robust

clean:
	rm -rf *.dat

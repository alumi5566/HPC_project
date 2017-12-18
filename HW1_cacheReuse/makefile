
all:
	gcc dgemm_3_2_1_0.c -o dgemm_3_2_1_0
	gcc dgemm_2_1_0.c -o dgemm_2_1_0
	gcc block_1_2_ijk_jik.c -o block_1_2_ijk_jik
	gcc block_1_2_jki_kji.c -o block_1_2_jki_kji
	gcc block_1_2_kij_ikj.c -o block_1_2_kij_ikj
	gcc block_both.c -o block_both
	gcc block_both.c -o block_both_O0 -O0
	gcc block_both.c -o block_both_O1 -O1
	gcc block_both.c -o block_both_O2 -O2
	gcc block_both.c -o block_both_O3 -O3

clean:
	rm -f dgemm_3_2_1_0
	rm -f dgemm_2_1_0
	rm -f block_1_2_ijk_jik
	rm -f block_1_2_jki_kji
	rm -f block_1_2_kij_ikj
	rm -f block_both
	rm -f block_both_O0
	rm -f block_both_O1
	rm -f block_both_O2
	rm -f block_both_O3 

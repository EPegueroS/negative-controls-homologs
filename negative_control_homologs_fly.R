negative_control_homologs_fly = function() {

	system("perl busca_homologos_en_grupos_shuffling.pl /home/epeguero/motif_find_meme/MEME_8_8_opt/MOSCA/BLAST_mosca /home/epeguero/motif_find_meme/MEME_8_8_opt/MAST /home/epeguero/motif_find_meme/MEME_8_8_opt/MOSCA ./datos_combinados_secuencias.txt /home/epeguero/motif_find_meme/mosca.ptt 64 > homolog_hits.txt")
	

}

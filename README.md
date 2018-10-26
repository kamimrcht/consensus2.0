test de détection de la période avec l'autocorrelation en fonction du taux d'erreur :
python generate_noise.py -e 0;  python3.6 autocorr.py error_seq_test.fa 
=> 0% d'erreur, 1e séquence a 3 répétitions et une période de 2804, 2e période sans pattern

python generate_noise.py -e 1;  python3.6 autocorr.py error_seq_test.fa
python generate_noise.py -e 10;  python3.6 autocorr.py error_seq_test.fa

etc
Construction de la séquence de reference
(a+c+b+c'+a'+b)*3
où a est l'adaptor
c la séquence,
b le bell
' sont les rc

#meme chose mais avec des quasi kmers
python generate_noise.py -e 1;  python3.6 autocorr_quasi_kmer.py error_seq_test.fa




update 26/10
python generate_noise2.py -e 10;  python3.6 autocorr_quasi_kmer.py error_seq_test.fa ; python validate_consensus.py   ### donne une idée de la qualité de acbc'a'b


python compute_consensus_from_period.py ; python validate_final_consensus.py # valide le consensus pour 1 read

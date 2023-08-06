# find stop codons in fasta file with headers and '>'
with open('Eukprot_db/kineto_free_NT_DM') as file:
    s = file.readlines()
    s_count = 0
    for i in s[1::2]:
        s_count += 1
        codons = [i[j:j+3] for j in range(0, len(i), 3)]
        stops = ['TAG', 'TAA', 'TGA']
        for stop in stops:
            if stop in codons:
                print(s_count, codons.index(stop)*3)

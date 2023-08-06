# prepare a fasta output from macse to Datamonkey - delete stop-codons in the ends and change ! to ?
input_file = 'Eukprot_db/kineto_free_NT'
output = 'Eukprot_db/kineto_free_NT_DM'
with open(input_file) as file:
    s = file.read()
    s = s.replace('!', '?')
    print(s)
with open(output, 'w') as file:
    file.write(s)
# nt = 'ATGC'
stop = ['TAG', 'TGA', 'TAA']
with open(output) as file:
    s = file.readlines()
    s = [s[i] if i % 2 == 0 else s[i] for i in range(len(s))]
    s = [s[i] if (i % 2 == 0 or (s[i][s[i].rfind('T'):s[i].rfind('T') + 3] not in stop))
         else s[i][:s[i].rfind('T')] + '---' + s[i][s[i].rfind('T') + 3:] for i
         in range(len(s))]
    s = ''.join(s)
    print(s)
with open(output, 'w') as file:
    file.write(s)

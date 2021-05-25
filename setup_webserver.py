import os
from netifaces import interfaces, ifaddresses, AF_INET

def ip4_addresses():
    ip_list = []
    for interface in interfaces():
        try:
            for link in ifaddresses(interface)[AF_INET]:
                ip_list.append(link['addr'])
        except Exception:
            pass
    return ip_list

print('setup IP address')
iplist = ip4_addresses()
f = open('/home/ag-rossi/projects/nanodorf/NGS_webserver/settings.py', 'r').readlines()
f[27] = f[27][:-2]+','+ str(iplist)[1:] + '\n'
f1 = open('/home/ag-rossi/projects/nanodorf/NGS_webserver/settings.py', 'w')
f1.writelines(f)
f1.close() 
  
print('create human reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_human.mmi ~/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
print('create mouse reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_mouse.mmi ~/ReferenceData/Mus_musculus.GRCm39.dna.toplevel.fa.gz')
print('create rat reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_rat.mmi ~/ReferenceData/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz')
print('create zebrafish reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_zebrafish.mmi ~/ReferenceData/Danio_rerio.GRCz11.dna.toplevel.fa.gz')

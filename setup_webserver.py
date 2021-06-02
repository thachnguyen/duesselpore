#Virtual box: to shrink VM image 
#mount -n -o remount,ro -t ext2 /dev/sda1 / 
#zerofree -v /dev/sda1
#reboot and run on your host: VBoxManage modifymedium disk "new_hdd.vdi" --compact
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

print('Downloading human reference genome')
os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P ~/ReferenceData/')
print('Downloading mouse reference genome')
os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz -P ~/ReferenceData/')
print('Downloading rat reference genome')
os.system('wget ftp://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz -P ~/ReferenceData/')
print('Downloading zebrafish reference genome')
os.system('wget ftp://ftp.ensembl.org/pub/release-104/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz -P ~/ReferenceData/')
print('Downloading C elegans reference genome')
os.system('wget ftp://ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz -P ~/ReferenceData/') 
  
print('creating human reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_human.mmi ~/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
print('creating mouse reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_mouse.mmi ~/ReferenceData/Mus_musculus.GRCm38.dna.toplevel.fa.gz')
print('creating rat reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_rat.mmi ~/ReferenceData/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz')
print('creating zebrafish reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_zebrafish.mmi ~/ReferenceData/Danio_rerio.GRCz11.dna.toplevel.fa.gz')
print('creating C elegans reference genome indexes, please wait')
os.system('minimap2 -t 4 -k14 -w5 -d ~/ReferenceData/reference_celegans.mmi ~/ReferenceData/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz')


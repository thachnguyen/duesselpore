from netifaces import interfaces, ifaddresses, AF_INET
import os, sys


def ip4_addresses():
    ip_list = []
    for interface in interfaces():
        try:
            for link in ifaddresses(interface)[AF_INET]:
                ip_list.append(link['addr'])
        except Exception:
            pass
    return ip_list

if __name__=="__main__":
	print('setup IP address')
	iplist = ip4_addresses()
	f = open('/home/ag-rossi/projects/duesselpore/NGS_webserver/settings.py', 'r').readlines()
	f[27] = f[27][:-2]+','+ str(iplist)[1:] + '\n'
	f1 = open('/home/ag-rossi/projects/duesselpore/NGS_webserver/settings.py', 'w')
	f1.writelines(f)
	f1.close()

	if sys.argv[1]=='light':
		print('Downloading human reference genome')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P ~/ReferenceData/')
		os.system('wget http://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz -P ~/ReferenceData/')
		#print('creating human reference genome indexes, please wait')
		#os.system('minimap2 -t 4 -k15 -w10 -d ~/ReferenceData/reference_human.mmi ~/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
		#os.system('echo %s|sudo -S %s'%(sudo_pass, install_cmd))
		#os.unlink('~/ReferenceData/*.fa.gz')

	if sys.argv[1]=='full':
		print('Downloading human reference genome')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P ~/ReferenceData/')
		os.system('wget http://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz -P ~/ReferenceData/')
		# print('creating human reference genome indexes, please wait')
		# os.system('minimap2 -t 4 -k15 -w10 -d ~/ReferenceData/reference_human.mmi ~/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
		# install_cmd ='apt install fastqc'
		# sudo_pass = '123456'
		# os.system('echo %s|sudo -S %s'%(sudo_pass, install_cmd))
		# os.unlink('~/ReferenceData/*.fa.gz')


		print('Downloading mouse reference genome')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -P ~/ReferenceData/')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz -P ~/ReferenceData/')

		print('Downloading rat reference genome')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz -P ~/ReferenceData/')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz -P ~/ReferenceData/')

		print('Downloading zebrafish reference genome')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz -P ~/ReferenceData/')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/danio_rerio/Danio_rerio.GRCz11.102.gtf.gz -P ~/ReferenceData/')

		print('Downloading C elegans reference genome')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz -P ~/ReferenceData/')
		os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.102.gtf.gz -P ~/ReferenceData/') 
 
		# print('creating mouse reference genome indexes, please wait')
		# os.system('minimap2 -t 4 -k15 -w10 -d ~/ReferenceData/reference_mouse.mmi ~/ReferenceData/Mus_musculus.GRCm39.dna.toplevel.fa.gz')
		# print('creating rat reference genome indexes, please wait')
		# os.system('minimap2 -t 4 -k15 -w10 -d ~/ReferenceData/reference_rat.mmi ~/ReferenceData/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz')
		# print('creating zebrafish reference genome indexes, please wait')
		# os.system('minimap2 -t 4 -k15 -w10 -d ~/ReferenceData/reference_zebrafish.mmi ~/ReferenceData/Danio_rerio.GRCz11.dna.toplevel.fa.gz')
		# print('creating C elegans reference genome indexes, please wait')
		# os.system('minimap2 -t 4 -k15 -w10 -d ~/ReferenceData/reference_celegans.mmi ~/ReferenceData/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz')
		# os.unlink('~/ReferenceData/*.fa.gz')

	print('Your webserver IP address is %s, please use http://%s:8000/duesselpore/ on your browser:'%(iplist[-1], iplist[-1]))
rule download_remote:
	output: 'databases/remotefiles/{f}'
	params: 
		url=lambda wildcards: config['remotefiles'][wildcards.f]['url'],
		md5=lambda wildcards: config['remotefiles'][wildcards.f]['md5']
	shell:
		'''
		curl -L {params.url} > {output}
		echo {params.md5} {output} | md5sum -c -
		'''


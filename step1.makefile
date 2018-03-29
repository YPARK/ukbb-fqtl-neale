all: ldblocks/EUR

ldblocks/EUR: ldblocks
	ln -s ./nygcresearch-ldetect-data-ac125e47bf7f/EUR $@

ldblocks: ac125e47bf7f.zip
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	unzip $< -d $@

ac125e47bf7f.zip:
	wget https://bitbucket.org/nygcresearch/ldetect-data/get/ac125e47bf7f.zip


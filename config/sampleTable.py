import yaml

with open("ecoli.yaml", "r") as stream:
    try:
        obj = yaml.safe_load(stream)
        meta = obj['meta'] 
        samples = list(meta.keys())
        firstKey = samples[0]
        fields = list(meta[firstKey].keys())
        print('\t'.join(['no', 'sample'] + fields))
        i = 1
        for sample in samples:
            if i != 1:
                print('\n', end="")
            print(str(i) + '\t' + sample, end="")
            i += 1
            for field in fields:
                print('\t' + str(meta[sample][field]), end="")
        print()
    except yaml.YAMLError as exc:
        print(exc)
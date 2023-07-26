def bacterialProteomeDomainAnalyzer(file, arg_type = None):
    import json 
    """sumary_line
    a rapid implementation of the bacterial domain
    analyzer, following the predictions of the bacterial
    domains from the interproscan. I implemented a mapped 
    dataframe approach to make it faster and iterable. it will 
    parse a nested to nested json from interpro for direct analysis
    and fecthing all the protein domains and the corresponding start
    and stop coordinates. Also you can make a direct ingestion to the
    database and it also provides a dataframe.
    Keyword arguments:
    argument -- file prediction by the interproscan
    Return: a systematic prediction of the domains in the sequences
    """
    if arg_type == "sequence":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        sequence = ''.join([i["sequence"] for i  in data["results"]])
        return sequence
    if arg_type == "interpro_normalize":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        interpro_normalize = pd.concat(list(map(lambda n: pd.DataFrame(n), \
                                [i["matches"] for i  in data["results"]])))
        return interpro_normalize
    if arg_type == "signature":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        signature = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['signature']
        return signature
    if arg_type == "location":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        location = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['locations']
        return location
    if arg_type == "evalue":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        evalue = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['evalue']
        return evalue
    if arg_type == "score":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        score = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['score']
        return score
    if arg_type == "modelac":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        modelac = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['model-ac']
        return modelac
    if arg_type == "scope":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        scope = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                               [i["matches"] for i  in data["results"]])))['scope']
        return scope
    if arg_type == "accession":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        accession = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                               [i["matches"] for i  in data["results"]])))['accession']
        return accession
    if arg_type == "name":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        name = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['name']
        return name
    if arg_type == "proteinClass":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        proteinClass = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['ProteinClass']
        return proteinClass
    if arg_type == "graftPoint":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        graftPoint = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['graftPoint']
        return graftPoint
    if arg_type == "goxRefs":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        goxRefs = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['goxRefs']
        return goxRefs
    if arg_type == "prediction_locations":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        prediction_locations = pd.DataFrame.from_dict(pd.concat(list(map(lambda n:\
                            pd.DataFrame(n),[i["matches"] for i  in data["results"]])))\
                                                        ["signature"].apply(lambda n: n.values()))
        return prediction_locations
    if arg_type == "getdomains":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        get_domains = list(list(filter(lambda n: n!=None and n!="",i)) for i in \
                        (list(list(filter(lambda n: not isinstance(n,dict),i)) \
                                for i in ([list(i) for i in pd.DataFrame.from_dict \
                                            (pd.concat(list(map(lambda n: pd.DataFrame \
                                                (n),[i["matches"] for i  in data["results"]]))) \
                              ["signature"].apply(lambda n: n.values()))["signature"].to_list()]))))
        return get_domains

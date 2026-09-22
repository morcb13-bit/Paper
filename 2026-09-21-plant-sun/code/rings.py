import json,b13_chain_units as U
rows,place,offs=U.build_stack()
json.dump([[list(c) for c in r] for r in place],open('R14.json','w'))
print([len(r) for r in place])

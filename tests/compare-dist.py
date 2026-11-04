import json,sys
expected=list(map(float,open(sys.argv[1], encoding='utf-8').read().splitlines()))
for path,want in zip(sys.argv[2:],expected):
    got=json.load(open(path, encoding='utf-8'))["value"][0]
    assert float(got) == want, (path,got,want)
assert len(sys.argv)-2 == len(expected)

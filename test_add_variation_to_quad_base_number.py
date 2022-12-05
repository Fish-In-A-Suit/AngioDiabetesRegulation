import itertools

length = 10
input = [n*4**x for n,x in itertools.product([1,2,3], range(1,length+1))]
print(len(input))

combinations = list(itertools.combinations(input, r=5))

print(len(combinations))
print(combinations)
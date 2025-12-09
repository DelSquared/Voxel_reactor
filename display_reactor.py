red   = "\033[31m██\033[0m"
blue  = "\033[94m██\033[0m"
white = "██"

print('# = ', white, ', - = ', blue, ', @ = ', red)
print()

with open("reactor.txt") as f:
  print(f.read().replace("#", white).replace("-", blue).replace("@", red))
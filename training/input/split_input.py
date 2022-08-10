inputFilesList_o=open("input.txt", "r")
inputFiles = inputFilesList_o.read().splitlines()
inputFilesList_o.close()

index=0
for i, name in enumerate(inputFiles):
   if i % 10 == 0:
      f=open("input_%i.txt" %index, "a")
      index=index+1
   f.write("%s \n" %name)
   if i % 10 == 9:
      f.close()

      

inputFilesList_o=open("input_2018_PFCands.txt", "r")
inputFiles = inputFilesList_o.read().splitlines()
inputFilesList_o.close()

index=0
for i, name in enumerate(inputFiles):
   if i % 50 == 0:
      f=open("input_PFCands_%i.txt" %index, "a")
      index=index+1
   f.write("%s \n" %name)
   if i % 10 == 49:
      f.close()

      

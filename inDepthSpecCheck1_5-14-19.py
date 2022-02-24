import os, re, datetime

p = re.compile('(?<=range=)[^//]*(?=/)')
p1 = re.compile('(?<=5\'pad=)\d*(?=/)')
p2 = re.compile('(?<=3\'pad=)\d*')
#sorting key

def myKey(s):
	if s[-3:] == 'psl':
		return 
	else:
		return 0

timeStamp = datetime.datetime.now().strftime('%Y_%m_%d_%H:%M:%S')
outFolder = os.getcwd() + '/BatchAnalysis'
RunFolder = outFolder + '/Run_' + timeStamp
ProcessedFolder = os.getcwd() + '/ProcessedFilteredPsl/'
if not os.path.exists(outFolder):os.makedirs(outFolder)
if not os.path.exists(ProcessedFolder):os.makedirs(ProcessedFolder)
if not os.path.exists(RunFolder):os.makedirs(RunFolder)
SpecOutDir = RunFolder + '/Specific/'
SharedOutDir = RunFolder + '/Shared/'
CheckOutDir = RunFolder + '/CheckAgain/'
if not os.path.exists(SpecOutDir):os.makedirs(SpecOutDir)
if not os.path.exists(SharedOutDir):os.makedirs(SharedOutDir)
if not os.path.exists(CheckOutDir):os.makedirs(CheckOutDir)

def checkFolders():
	Spec, Shared, Check = [], [], []
	SpecFiles = os.listdir(SpecOutDir)
	pattern = re.compile('(?<=Specific)\d*(?=.psl)')
	SpecFiles = pattern.findall(str(SpecFiles))
	if len(SpecFiles) == 0:
		Spec = 0
	else:
		for num in SpecFiles:
			Spec.append(int(num))
		Spec = max(Spec) + 1
	SharedFiles = os.listdir(SharedOutDir)
	pattern = re.compile('(?<=Shared)\d*(?=.psl)')
	SharedFiles = pattern.findall(str(SharedFiles))
	if len(SharedFiles) == 0:
		Shared = 0
	else:
		for num in SharedFiles:
			Shared.append(int(num))
		Shared = max(Shared) + 1
	CheckFiles = os.listdir(CheckOutDir)
	pattern = re.compile('(?<=CheckAgain)\d*(?=.psl)')
	CheckFiles = pattern.findall(str(CheckFiles))
	if len(CheckFiles) == 0:
		Check = 0
	else:
		for num in CheckFiles:
			Check.append(int(num))
		Check = max(Check) + 1
	return Spec, Shared, Check

SpecNumber, SharedNumber, CheckAgainNum = checkFolders()

indir = os.getcwd() + '/BatchPslFiltered/'

	# PSL COLUMNS --> [0] = Match, [1] = mismatch, [2] = rep. Match, [3] = N's, [4] = Q Gap count, [5] = Qgap bases
	# [6] = T Gap count, [7] = T gap bases, [8] = Strand, [9] = Q name, [10] = Q size, [11] = Q Start
	#[12] = Q end, [13] = T name, [14] = Tsize, [15] = T start, [16] = T end, [17] = BlockCount, [18] = Block Sizes
	# [19] = qStarts, [20] = tStarts

filenames = os.listdir(indir)
print(filenames)
total = 0
preSpecNumber = SpecNumber
preSharedNumber = SharedNumber
preCheckAgainNum = CheckAgainNum
for f in sorted(filenames, key=myKey):
	psl = open(os.path.join(indir, f), 'r')
	for line in psl.readlines():
		total += 1
		columns = line.split('\t')
		leftPad = int(p1.findall(columns[9])[0])
		rightPad = int(p2.findall(columns[9])[0])
		qStarts = columns[19].split(',')
		blockSizes = columns[18].split(',')
		tStarts = columns[20].split(',')
		qGapStarts = []
		qGapSizes = []
		qGapDistanceFromStart = []
		tGapStarts = []
		tGapSizes = []
		del qStarts[-1]
		del blockSizes[-1]
		del tStarts[-1]

		for n in range(len(qStarts)-1):
			Gap = int(qStarts[n]) + int(blockSizes[n])
			Size = int(qStarts[n+1]) - Gap
			qGapStarts.append(Gap)
			qGapSizes.append(Size)
			qGapDistanceFromStart.append(abs(int(qGapStarts[-1]) - leftPad))
			
		for n in range(len(tStarts)-1):
			Gap = int(tStarts[n]) + int(blockSizes[n])
			Size = int(tStarts[n+1]) - Gap
			tGapStarts.append(Gap)
			tGapSizes.append(Size)
		try:
			index = qGapDistanceFromStart.index(min(qGapDistanceFromStart))
			aluStart = qGapStarts[index]
			aluSize = qGapSizes[index]
		except ValueError as e:
			aluStart = 1
			aluSize = 1

		if leftPad - 60 < aluStart <= rightPad + 20 and aluSize > 230:
			if tGapSizes[index] < 50:
				outName = SpecOutDir + '/Specific' + str(SpecNumber) + '.psl'
				out = open(outName, 'w')
				out.write(line)
				out.close
				SpecNumber = SpecNumber + 1
			elif abs(tGapSizes[index] - qGapSizes[index]) < 20 or tGapSizes[index] > 350:
				outName = SharedOutDir + '/Shared' + str(SharedNumber) + '.psl'
				out = open(outName, 'w')
				out.write(line)
				out.close
				SharedNumber = SharedNumber + 1
			elif tGapSizes.index(max(tGapSizes)) != index and max(tGapSizes) > 200:
				outName = CheckOutDir + '/CheckAgain' + str(CheckAgainNum) + '.psl'
				out = open(outName, 'w')
				out.write(line)
				out.close
				CheckAgainNum = CheckAgainNum + 1
			else:
				outName = CheckOutDir + '/CheckAgain' + str(CheckAgainNum) + '.psl'
				out = open(outName, 'w')
				out.write(line)
				out.close
				CheckAgainNum = CheckAgainNum + 1
		else:
			outName = SharedOutDir + '/Shared' + str(SharedNumber) + '.psl'
			out = open(outName, 'w')
			out.write(columns[9] + ' ---->\tN\'s - [' + columns[3] + ']\tQGapBases - [' + columns[5] + ']\tTargetGapBases - [' + columns[7] + ']\n\n' + line)
			out.close
			SharedNumber = SharedNumber + 1
	os.rename(psl.name, ProcessedFolder + psl.name.split('/')[-1])

# lines = ['\n# of sequences analyzed:  ', str(total), '\n# of specific sequences:  ', \
# str(SpecNumber - preSpecNumber), ' (', str((SpecNumber - preSpecNumber) / total*100), '%)', \
# '\n# of shared sequences:  ', str(SharedNumber - preSharedNumber), ' (', str((SharedNumber - preSharedNumber) / total*100), '%)', \
# '\n# of sequences that were inconclusive:  ', str(CheckAgainNum - preCheckAgainNum), ' (', str((CheckAgainNum - preCheckAgainNum) / total*100), '%)']
# logDir = '/home/thomas/pythonFiles/primerDesign/logFiles/'
# logFiles = os.listdir(logDir)
# log = open(os.path.join(logDir, sorted(logFiles, reverse=True)[0]), 'a')
# log.writelines(lines)
# log.close()


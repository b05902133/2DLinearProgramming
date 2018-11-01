# PROG from ../makefile
program := ../${PROG}

testcaseDir   := Testcase
testResultDir := TestResult
testcases     := $(notdir $(wildcard ${testcaseDir}/*.input))
testResults   := $(addprefix ${testResultDir}/,$(subst .input,.log,${testcases}))
testRefs      := $(addprefix ${testcaseDir}/,$(subst .input,.output,${testcases}))

timeLimit := 1     # s
memLimit  := 64000 # KB

.PHONY: all clean ${testRefs}

all: ${testResultDir} ${testResults} ${testRefs}

${testResultDir}:
	mkdir $@

${testRefs}: testResult += $(patsubst ${testcaseDir}/%.output,${testResultDir}/%.log,$@)
${testRefs}:
	@diff $@ ${testResult} > /dev/null; \
	if [ $$? -eq 0 ]; then \
	  echo "${testResult} test pass!"; \
	else \
	  echo "${testResult} test fail!"; \
	fi

${testResults}: testcase += $(patsubst ${testResultDir}/%.log,${testcaseDir}/%.input,$@)
${testResults}: ${program}
	ulimit -t ${timeLimit} -v ${memLimit}; \
	${program} > $@ < ${testcase}

clean:
	rm -r ${testResultDir}

#!/usr/bin/env ruby

inputFileName = ARGV[0]
outputFileName = inputFileName + "result"
puts inputFileName
puts outputFileName

output = File.open(outputFileName,'w')

File.open(inputFileName, 'r') do |file|
  file.each_line do |line|
    tokens = line.split()
    if tokens.length == 3
      if tokens[0] == "PreSignRate"
        output.print tokens[2]
      end
      if tokens[0] == "PostSignRate"
        output.print ' ', tokens[2], "\n"
      end
    end
  end
end

output = File.open('174.data','w')

File.open("174.log", 'r') do |file|
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

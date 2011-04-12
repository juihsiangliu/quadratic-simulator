NAME = 'vin_10.txt';
V = 1.2;
STEP = 10;

% ==================================
inc = V/(STEP/2);
vin = V*ones(STEP,1);
j = 1;
for(i=0:inc:V)
	vin(j) = i;
	j = j+1;
end
% ==================================
fp = fopen(NAME,'w');
for(i=1:STEP-1)
	fprintf(fp,'%g,',vin(i));
end
fprintf(fp,'%g',vin(STEP));
fclose(fp);



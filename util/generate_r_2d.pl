$size = 100;
$nodeNum = $size*$size;
$endTime = "0.000000000001";
$timeStep = "0.000000000001";

$vinM = "vin_pwl_m.txt";
$vinAlpha = "vin_zero_1500_alpha.txt";
$vinBeta = "vin_zero_1500_beta.txt";
$vinGamma = "vin_zero_1500_gamma.txt";

$gAlpha = "g_alpha.txt";
$gBeta = "g_beta.txt";
$gGamma = "g_gamma.txt";

$cAlpha = "c_alpha.txt";
$cBeta = "c_beta.txt";
$cGamma = "c_gamma.txt";

$gVal = "100";

printf("*global_var_num\n");
printf("*2\n");
printf("*global_cor_matrix\n");
printf("*1.0 0.0\n");
printf("*0.0 1.0\n");
printf("*node_num\n");
printf("*%d\n",$nodeNum);
printf("*end_time\n");
printf("*%s\n",$endTime);
printf("*time_step\n");
printf("*%s\n",$timeStep);

printf("V10 1 0 vin_pwl_m.txt vin_zero_1500_alpha.txt vin_zero_1500_beta.txt vin_zero_1500_gamma.txt\n");

$index = 1;
for($i=1;$i<=$size;$i=$i+1)
{
	for($j=1;$j<$size;$j=$j+1)
	{
		$pos = ($i-1)*$size + $j;
		$neg = $pos+1;
		printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
		$index = $index+1;
	}
}


for($i=1;$i<$size;$i=$i+1)
{
	for($j=1;$j<=$size;$j=$j+1)
	{
		$pos = ($i-1)*$size + $j;
		$neg = $pos + $size;
		printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
		$index = $index+1;
	}
}



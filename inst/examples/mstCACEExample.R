if(interactive()){
  
data(mstData)

######################## weighted ITT ##############################################

caceOutput<- caceMSTBoot(Posttest~ Prettest+ Intervention,
			random="School",intervention="Intervention",
			compliance = "Percentage_Attendance",nBoot=1000,data=mstData)

cace <- caceOutput$CACE
cace

Complier <- caceOutput$Compliers
Complier 

### visualising CACE effect size

plot(caceOutput)
}

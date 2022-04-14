#!/bin/sh

local r
r=0

for i in CARDAH CARDAM CARDA_OM CARDA-OM CARDAW CARDG CARDOPTO E_dhtom E_DYN_EQ EXTRACT FRAME4P FRAME4V FRAMEP FRAMEV GTOM INTERSLP JTOJ MAKEL0 MAKEL MAKEL_P MTOSCREW PROJPONL ROTAT RTOCARDA SCREWTOM TRASF_MA TRASF_MH TR_MAMT VELWH2 WTOL_P WTOL_R WTOVEL WTOVEL_P ; do
	LD_LIBRARY_PATH="$LD_LIBRARY_PATH:../../" ./"${i}" > "${i}.temp.txt"
	if cmp -s "${i}.temp.txt" "${i}.txt"; then
		echo -e "SUCCESSO:\t${i}"
	else
		echo -e "FALLITO:\t${i}\nOTTENUTO:"
		echo "------------------------------------------------------------------------"
		cat "${i}.temp.txt"
		echo "------------------------------------------------------------------------"
		echo -e "INVECE DI:"
		echo "------------------------------------------------------------------------"
		cat "${i}.txt"
		echo "------------------------------------------------------------------------"
		echo -e "\n"
		r=-1
	fi
	rm "${i}.temp.txt"
done

return ${r}

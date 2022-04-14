#!/bin/sh

r=0

for i in CARDAH CARDAM CARDA_OM CARDA-OM CARDAW CARDG CARDOPTO E_dhtom E_DYN_EQ EXTRACT FRAME4P FRAME4V FRAMEP FRAMEV GTOM INTERSLP JTOJ MAKEL0 MAKEL MAKEL_P MTOSCREW PROJPONL ROTAT RTOCARDA SCREWTOM TRASF_MA TRASF_MH TR_MAMT VELWH2 WTOL_P WTOL_R WTOVEL WTOVEL_P ; do
	LD_LIBRARY_PATH="$LD_LIBRARY_PATH:../../" ./"${i}" > "${i}.tmp"

	echo "------------------------------------------------------------------------"
	cat "${i}.tmp"
	echo "------------------------------------------------------------------------"
	echo "------------------------------------------------------------------------"
	cat "${i}.txt"
	echo "------------------------------------------------------------------------"

	if cmp -s "${i}.tmp" "${i}.txt"; then
		echo -e "SUCCESSO:\t${i}"
	else
		echo -e "FALLITO:\t${i}"
		r=-1
	fi
done

exit ${r}

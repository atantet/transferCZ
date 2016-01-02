load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

pi = 3.14159265
REarth = 6378137
m2deg = 180. / (REarth * pi)
  
begin

  dataDir = "$HOME/PhD/data/"
  dstGridDir = "$HOME/PhD/dev/ENSO/init/"
  dstDir = "./"
  project = "FSU"
  dataset = "PAC"
  srcFile = "pac_pstress_6199.nc"
  srcVar1 = "Wu"
  srcVar2 = "Wv"
  dstNLat = 31
  dstNLon = 30
  startLon = 141.
  year0 = 1961
  yearEnd = 1994
  startt = 0
  
  srcDir = dataDir + project + "/" + dataset + "/"
  nYears = yearEnd - year0 + 1
  nt = nYears * 12

                                ; Read destination coordinates
  dstGridName = sprinti("%d", dstNLat) + "x" + sprinti("%d", dstNLon)
  dstLat = asciiread(dstGridDir + "y_" + dstGridName + ".txt", -1, "float") * m2deg
  dstLon = asciiread(dstGridDir + "x_" + dstGridName + ".txt", -1, "float") * m2deg + startLon
  asciiwrite(dstGridDir + "lat_" + dstGridName + ".txt", dstLat)
  asciiwrite(dstGridDir + "lon_" + dstGridName + ".txt", dstLon)
  
; Read first netcdf file to get coordinates
  srcPath = srcDir + srcFile
  dset = addfile(srcPath, "r")
  srcData1 = dset->$srcVar1$(startt:startt+nt-1, :, :)
  srcData2 = dset->$srcVar2$(startt:startt+nt-1, :, :)
  srcLat = dset->lat
  srcLon = dset->lon

  ; Allocate
  dstData1 = new((/nt, dstNLat, dstNLon/), "float")
  dstData2 = new((/nt, dstNLat, dstNLon/), "float")
  time = fspan(0, nt - 1, nt)

  do t = 0, nt - 1
                                ; interpolate
    dstData1(t, :, :) = linint2(srcLon, srcLat, srcData1(t, :, :), True, dstLon, dstLat, 0) * srcData1@scale_factor
    dstData2(t, :, :) = linint2(srcLon, srcLat, srcData2(t, :, :), True, dstLon, dstLat, 0) * srcData2@scale_factor
  end do

  dstData1!0 = "time"      
  dstData1!1 = "lat"      
  dstData1!2 = "lon"      
  dstData1&time = time
  dstData1&lat = dstLat
  dstData1&lon = dstLon
  dstData2!0 = "time"      
  dstData2!1 = "lat"      
  dstData2!2 = "lon"      
  dstData2&time = time
  dstData2&lat = dstLat
  dstData2&lon = dstLon

  dstFile = "pac." + dstGridName + "_" + sprinti("%d", year0) + sprinti("%d", yearEnd) + ".nc"
  print("Writing " + dstDir + dstFile)
  if (isfilepresent(dstDir + dstFile)) then
    system("rm " + dstDir + dstFile)
  end if
  dset = addfile(dstDir + dstFile, "c")
  dset->$srcVar1$ = dstData1
  dset->$srcVar2$ = dstData2
  delete(dset)
end
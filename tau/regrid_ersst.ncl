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
  project = "NCDC"
  dataset = "ERSST"
  prefix = "ersst."
  srcVar = "sst"
  dstNLat = 31
  dstNLon = 30
  startLon = 141.
  year0 = 1961
  yearEnd = 1994
  
  srcDir = dataDir + project + "/" + dataset + "/"
  nYears = yearEnd - year0 + 1
  nt = nYears * 12

  ; Read destination coordinates
  dstGridName = sprinti("%d", dstNLat) + "x" + sprinti("%d", dstNLon)
  dstLat = asciiread(dstGridDir + "lat_" + dstGridName + ".txt", -1, "float") * m2deg
  dstLon = asciiread(dstGridDir + "lon_" + dstGridName + ".txt", -1, "float") * m2deg + startLon
  
; Read first netcdf file to get coordinates
  srcFile = srcDir + prefix + sprinti("%d", year0) + "01.nc"
  dset = addfile(srcFile, "r")
  srcLat = dset->lat
  srcLon = dset->lon
  delete(dset)

  ; Allocate
  dstData = new((/nt, dstNLat, dstNLon/), "float")
  time = fspan(0, nt - 1, nt)

  do y = 0, nYears - 1
    year = year0 + y
    do m = 0, 11
      month = m + 1
      ; Read record
      srcFile = srcDir + prefix + sprinti("%d", year) + sprinti("%02d", month) + ".nc"
      print("Reading " + srcFile)
      dset = addfile(srcFile, "r")
      srcData = dset->$srcVar$(0, 0, :, :)
      ; interpolate
      dstData(y*12+m, :, :) = linint2(srcLon, srcLat, srcData, True, dstLon, dstLat, 0) * srcData@scale_factor + srcData@add_offset
      delete(srcData)
      delete(dset)
    end do
  end do

  dstData!0 = "time"      
  dstData!1 = "lat"      
  dstData!2 = "lon"      
  dstData&time = time
  dstData&lat = dstLat
  dstData&lon = dstLon

  dstFile = prefix + dstGridName + "_" + sprinti("%d", year0) + sprinti("%d", yearEnd) + ".nc"
  print("Writing " + dstDir + dstFile)
  if (isfilepresent(dstDir + dstFile)) then
    system("rm " + dstDir + dstFile)
  end if
  dset = addfile(dstDir + dstFile, "c")
  dset->$srcVar$ = dstData
  delete(dset)

end
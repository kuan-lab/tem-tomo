#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# make_nloo_joblist.sh  – generate a self-contained bash script (nloo_jobs.sh)
#                        with one et-nloo command per sub-volume.
# ---------------------------------------------------------------------------
set -euo pipefail

### ---------- USER PARAMETERS ------------------------------------------------
num_cores=8                 # only used for the helper message
cube_size=45

region_x=1891               # centred window sizes  ( -1 ⇒ full axis )
region_y=1891
region_z=136

skip_x=4                  # skip factors (cube multiples)
skip_y=4
skip_z=1

imod_cmd="tilt_21_lim30.com"  
et_nloo_bin="$HOME/Tomo/electra_0.5.4/electra/bin/et-nloo"

out_dir="nloo2d_subvols_21_lim30"
mkdir -p "${out_dir}"
jobfile="${out_dir}/nloo_jobs.sh"
### --------------------------------------------------------------------------

# ---------- find reconstruction .mrc from OutputFile line -------------------
mrc_file=$(grep -i -m1 '^ *OutputFile .*\.mrc' "${imod_cmd}" \
           | awk '{for(i=1;i<=NF;i++) if($i~/\.mrc$/){print $i; break}}')
[[ -z $mrc_file ]] && { echo "ERROR: no OutputFile *.mrc found"; exit 1; }

# ---------- tomogram size ---------------------------------------------------
read NY NZ NX <<<"$(header -size "${mrc_file}")"

# ---------- centre window bounds -------------------------------------------
[[ $region_x -lt 0 ]] && region_x=$NX
[[ $region_y -lt 0 ]] && region_y=$NY
[[ $region_z -lt 0 ]] && region_z=$NZ

cx=$((NX/2)); cy=$((NY/2)); cz=$((NZ/2))
min_x=$(( cx - region_x/2 )); max_x=$(( cx + region_x/2 - cube_size ))
min_y=$(( cy - region_y/2 )); max_y=$(( cy + region_y/2 - cube_size ))
min_z=$(( cz - region_z/2 )); max_z=$(( cz + region_z/2 - cube_size ))

# clip to volume
min_x=$(( min_x < 0 ? 0 : min_x ));    max_x=$(( max_x > NX-cube_size ? NX-cube_size : max_x ))
min_y=$(( min_y < 0 ? 0 : min_y ));    max_y=$(( max_y > NY-cube_size ? NY-cube_size : max_y ))
min_z=$(( min_z < 0 ? 0 : min_z ));    max_z=$(( max_z > NZ-cube_size ? NZ-cube_size : max_z ))

# ---------- step sizes ------------------------------------------------------
step_x=$(( cube_size * skip_x ))
step_y=$(( cube_size * skip_y ))
step_z=$(( cube_size * skip_z ))

# ---------- how many cubes? -------------------------------------------------
nx=$(( (max_x - min_x) / step_x + 1 ))
ny=$(( (max_y - min_y) / step_y + 1 ))
nz=$(( (max_z - min_z) / step_z + 1 ))
(( max_x < min_x )) && nx=0
(( max_y < min_y )) && ny=0
(( max_z < min_z )) && nz=0
total=$(( nx * ny * nz ))

echo "=== NLOO grid summary ======================================="
echo "Tomogram size: NX=${NX}  NY=${NY}  NZ=${NZ}"
echo "Cube size:     ${cube_size}"
echo "Region ranges: X:[${min_x}-${max_x}]  Y:[${min_y}-${max_y}]  Z:[${min_z}-${max_z}]"
echo "Step sizes:    X=${step_x}  Y=${step_y}  Z=${step_z}"
echo "Cubes along:   nx=${nx}  ny=${ny}  nz=${nz}"
echo "TOTAL cubes:   ${total}"
echo "Job list will be written to ${jobfile}"
echo "-------------------------------------------------------------"

[[ $total -eq 0 ]] && { echo "ERROR: zero cubes fit – adjust parameters."; exit 1; }

# ---------- write job list --------------------------------------------------
printf '#!/usr/bin/env bash\n' > "$jobfile"
count=0
for (( z=min_z; z<=max_z; z+=step_z )); do
  for (( y=min_y; y<=max_y; y+=step_y )); do
    for (( x=min_x; x<=max_x; x+=step_x )); do
      (( ++count ))
      out_base="nloo2d_xyz_${x}_${y}_${z}"
      in_txt="${out_base}_nloo2d_res.txt"
      out_txt="${out_base}_nloo2d_res_est.txt"

      printf '%s\n' \
      "echo \"[#$(printf '%04d' "$count")] cube @ (${x},${y},${z})\" >&2 && \
       ${et_nloo_bin} -box ${x},${y},${z},${cube_size},${cube_size},${cube_size} \
               -imod ${imod_cmd} -nloo2d ${out_base} && \
       ${et_nloo_bin/et-nloo/et-tiltnloo} ${in_txt} ${out_txt} && \
       mv ${in_txt}  ${out_dir}/${in_txt}  && \
       mv ${out_txt} ${out_dir}/${out_txt}" \
       >> "$jobfile"
      
      
      #printf 'echo "[#%04d] cube @ (%d,%d,%d)" >&2\n'  "$count" "$x" "$y" "$z"  >> "$jobfile"
      #printf '%q -box %d,%d,%d,%d,%d,%d -imod %q -nloo2d %q\n' \
      #       "$et_nloo_bin" "$x" "$y" "$z" "$cube_size" "$cube_size" "$cube_size" \
      #       "$imod_cmd"    "$out_base" >> "$jobfile"
      
      ### --- SECOND job list: convert to tilt-resolution -----------------
      #conv_file="convert_jobs.sh"
      #if (( count == 1 )); then                 # first cube? ⇒ create & she-bang
      #    printf '#!/usr/bin/env bash\n' > "${out_dir}/$conv_file"
      #fi

      #in_txt="${out_base}_nloo2d_res.txt"
      #out_txt="${out_base}_nloo2d_res_est.txt"

      #printf '%s\n' \
      #"echo \"[conv #$(printf '%04d' "$count")] ${in_txt}\" >&2" \
      #"${et_nloo_bin/et-nloo/et-tiltnloo} ${in_txt} ${out_txt}" \
      #>> "$jobfile"

      # move output files into output subdir (nloo default writes to main dir)
      #printf '%s\n' "mv $in_txt $out_dir/$in_txt" >> "$jobfile"
      #printf '%s\n' "mv $out_txt $out_dir/$out_txt" >> "$jobfile" 

    done
  done
done

chmod +x "$jobfile" 
echo "Wrote $count et-nloo commands to $jobfile"
echo
echo "cd $out_dir"
echo "Run the recon jobs :  parallel -j $num_cores --line-buffer < $jobfile"

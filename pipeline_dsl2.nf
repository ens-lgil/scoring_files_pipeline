#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def list_to_string(data_list) {
  return data_list.join(',')
}


process get_variants_list {
  input:
    val pgs_ids

  output:
    val params.var_file_name

  script:
  """
  python $params.loc_pipeline/variants_list.py \
    --scores_ids $pgs_ids \
    --scores_dir $params.work_dir/scorefiles/ \
    --var_file $params.work_dir/variants/${params.var_file_name}.txt \
    --coord_file $params.work_dir/locations/$params.coord_file_name
  """
}


process var2location {
  input:
    val file_name

  output:
    val params.new_coord_file_name

  script:
  """
  perl $params.loc_harmonizer/EnsemblMappings/var2location_3738.pl $file_name $params.work_dir
  cp $params.work_dir/${file_name}.out $params.work_dir/locations/${params.new_coord_file_name}
  """
}


process update_varlocation_file {
  input:
    val new_file_name

  script:
  """
  python $params.loc_pipeline/update_variant_location_file.py \
    --loc_file $params.work_dir/locations/${new_file_name} \
    --coord_file $params.work_dir/locations/${params.coord_file_name}
  """
}


process validate_scoring_file {
  input:
    val pgs_id

  output:
    val pgs_id

  script:
  """
  rm -f $params.work_dir/logs/${pgs_id}_log.txt \
  python $params.loc_validator/run_validator.py \
    -f $params.work_dir/scorefiles/${pgs_id}.txt.gz \
    --log_dir $params.work_dir/logs/
  """
}


process validation_result {
  input:
    val pgs_id

  output:
    stdout emit: result_count

  script:
  """
  grep -m 1 "File is valid" $params.work_dir/logs/${pgs_id}_log.txt | wc -l
  """
}


process HmPOS {
  input:
    val pgs_id
    val is_valid

  output:
    val pgs_id

  when:
    is_valid =~ /1/

  script:
  """
  python $params.loc_harmonizer/Harmonize.py HmPOS $pgs_id $params.genebuild -loc_files $params.loc_score_dir -loc_hmoutput $params.loc_hmoutput --gzip
  """
}


process HmVCF {
  input:
    val pgs_id

  script:
  """
  python $params.loc_harmonizer/Harmonize.py HmVCF $pgs_id $params.genebuild -loc_files $params.loc_hmoutput -loc_hmoutput $params.loc_hmoutput -loc_vcfs $params.loc_vcfs
  """
}


workflow {
  data = channel.from(params.pgs)

  // Prepare variants list and their locations
  get_variants_list(list_to_string(params.pgs))
  var2location(get_variants_list.out)
  update_varlocation_file(var2location.out)

  // Validate scoring files
  validate_scoring_file(data)
  validation_result(validate_scoring_file.out)

  // Harmonise scoring files
  HmPOS(validate_scoring_file.out, validation_result.out)
  HmVCF(HmPOS.out)
}
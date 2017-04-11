# staircase
rm data_sfari_batch1[ABCD]*/*/*smlcsm*

# collapse
rm data_sfari_batch1[ABCD]*/*/*smlx2csm*

# priors
for i in data_sfari_batch1[ABCD]*; do
  psql -a -c "drop table ${i}.prior cascade; drop table ${i}.prior_region cascade"
done

# posterior
rm data_sfari_batch1[ABCD]*/*/*bayescsm*

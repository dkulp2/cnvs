# staircase
rm data_sfari_batch1[ABCD]*/*/*smlcsm*

# collapse
rm data_sfari_batch1[ABCD]*/*/*smlx2csm*

# priors and posterior tables - SET ENVIRONMENT!!
for i in data_sfari_batch1[ABCD]*; do
    psql -c "drop table ${i}.prior cascade; drop table ${i}.prior_region cascade"
    psql -c "drop table ${i}.posterior cascade; drop table ${i}.posterior_dist cascade;"
done

# posterior
rm data_sfari_batch1[ABCD]*/*/*bayescsm*

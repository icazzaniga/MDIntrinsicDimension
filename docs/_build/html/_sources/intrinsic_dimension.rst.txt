intrinsic_dimension()
============================

Estimates the intrinsic dimension (ID) of the object over an entire MD system.



.. code-block:: python

   from intrinsic_dimension import intrinsic_dimension

   #estimate local ID of villin
   mean_all, mean_last, local_id = intrinsic_dimension(topology='villin/2F4K.pdb', trajectory='villin/2F4K_F.xtc', projection_method='Dihedrals', id_method='TwoNN')

   #estimate global ID of villin
   global_all, global_last = intrinsic_dimension(topology='villin/2F4K.pdb', trajectory='villin/2F4K_F.xtc', projection_method='Dihedrals', id_method='TwoNN')

   #results local
   print('Mean instantaneous ID of the entire trajectory:', mean_all)
   print('Mean instantaneous ID of the last 100 frames:', mean_last)
   print('Istantaneous ID of the entire trajectory:', local_id)
   #results global
   print('Average ID of the entire trajectory:', global_all)
   print('Average ID of the last 100 frames:', global_last)

.. note::
   INSERIRE OUTPUT VILLINA
   ESEMPIO PLOT?


  










%parameter_hash = (
                    'param' => [
                                 'theta',
                                 'theta',
                                 'theta',
                                 'theta',
                                 'theta',
                                 'theta',
                                 'omega',
                                 'omega',
                                 'sigma'
                               ],
                    'sdcorrform' => [
                                      0,
                                      0,
                                      0,
                                      0,
                                      0,
                                      0,
                                      0,
                                      0,
                                      0
                                    ],
                    'upper_bounds' => [
                                        3,
                                        7,
                                        '0.8',
                                        5,
                                        8,
                                        '-0.001',
                                        1000000,
                                        1000000,
                                        1000000
                                      ],
                    'block_number' => [
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0
                                      ],
                    'inits' => [
                                 '1.41',
                                 '4.73',
                                 '0.231',
                                 '2.15',
                                 '5.78',
                                 '-0.127',
                                 '0.608',
                                 '0.375',
                                 '0.289'
                               ],
                    'coords' => [
                                  '1',
                                  '2',
                                  '3',
                                  '4',
                                  '5',
                                  '6',
                                  '1,1',
                                  '2,2',
                                  '1,1'
                                ],
                    'labels' => [
                                  '[1]   CL (L/d)',
                                  '[2,3] V1  (L)',
                                  'V2',
                                  'Q',
                                  'CLADA',
                                  'V1ALB1',
                                  'IIV CL',
                                  'IIV V1',
                                  'EPS(1), proportional'
                                ],
                    'lower_bounds' => [
                                        '0.001',
                                        '0.001',
                                        '0.01',
                                        '0.001',
                                        '0.001',
                                        '-0.140',
                                        0,
                                        0,
                                        0
                                      ],
                    'values' => [
                                  '1.39726',
                                  '4.70664',
                                  '0.228756',
                                  '2.1274',
                                  '5.78029',
                                  '-0.127136',
                                  '0.607047',
                                  '0.375377',
                                  '0.289024'
                                ],
                    'off_diagonal' => [
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0
                                      ],
                    'coordinate_strings' => [
                                              'THETA1',
                                              'THETA2',
                                              'THETA3',
                                              'THETA4',
                                              'THETA5',
                                              'THETA6',
                                              'OMEGA(1,1)',
                                              'OMEGA(2,2)',
                                              'SIGMA(1,1)'
                                            ],
                    'choleskyform' => [
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0
                                      ]
                  );
@center_rawresults_vector = (
                              'center',
                              25,
                              1,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              undef,
                              undef,
                              undef,
                              0,
                              undef,
                              undef,
                              'FOCE_I',
                              undef,
                              0,
                              0,
                              '10624.102283304297',
                              '1.40369',
                              '4.80049',
                              '0.175317',
                              '1.59169',
                              '5.69722',
                              '-0.12651',
                              '0.652189',
                              '0.3732',
                              '0.290549',
                              '1e-13',
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef,
                              undef
                            );
@negative_dofv = (
                   0,
                   0,
                   0,
                   0,
                   2
                 );
@intermediate_raw_results_files = (
                                    'raw_results_sir_iteration1.csv',
                                    'raw_results_sir_iteration2.csv',
                                    'raw_results_sir_iteration3.csv',
                                    'raw_results_sir_iteration4.csv',
                                    'raw_results_run12.csv'
                                  );
@samples = (
             1000,
             1000,
             1000,
             2000,
             2000
           );
@resamples = (
               200,
               400,
               500,
               1000,
               1000
             );
@attempted_samples = (
                       1000,
                       1000,
                       1000,
                       2000,
                       2000
                     );
@successful_samples = (
                        1000,
                        1000,
                        1000,
                        2000,
                        2000
                      );
@actual_resamples = (
                      200,
                      400,
                      500,
                      1000,
                      1000
                    );
@seed_array = (
                2196596,
                0
              );
$with_replacement = 0;
$mceta = 0;
$problems_per_file = 100;
@minimum_ofv = (
                 '10624.323128372396',
                 '10624.375740605170',
                 '10624.473486991717',
                 '10624.341891115700',
                 '10624.102283304297'
               );
$reference_ofv = '10624.155290760960';
$done = 1;
$iteration = 5;
$recenter = 1;
$copy_data = 1;
$boxcox = 1;
$nm_version = 'default';
$model_filename = 'run12.mod';
$cap_resampling = 1;
$cap_correlation = '0.8';
$adjust_blocks = 0;
$check_cholesky_reparameterization = 0;
$subjects = '75';

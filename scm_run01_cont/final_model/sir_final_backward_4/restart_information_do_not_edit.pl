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
                    'labels' => [
                                  'CL (L/d) [1]',
                                  'V (L) [2,3]',
                                  'CLADA',
                                  'VADA',
                                  'VALB',
                                  'VCREAT',
                                  'IIV CL (value obtained from SCM+)',
                                  'IIV V',
                                  'EPS(1), proportional (value obtained from SCM+)'
                                ],
                    'lower_bounds' => [
                                        '0.001',
                                        '0.001',
                                        '0.001',
                                        -1,
                                        '-0.138',
                                        '-0.4',
                                        0,
                                        0,
                                        0
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
                                 '1.4',
                                 '4.96',
                                 '4.02',
                                 '-0.263',
                                 '-0.125',
                                 '0.184',
                                 '0.604',
                                 '0.37',
                                 '0.283'
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
                                      ],
                    'upper_bounds' => [
                                        8,
                                        15,
                                        9,
                                        2,
                                        '0.084',
                                        '0.4',
                                        1000000,
                                        1000000,
                                        1000000
                                      ],
                    'values' => [
                                  '1.39797',
                                  '4.96413',
                                  '4.01979',
                                  '-0.262909',
                                  '-0.124656',
                                  '0.183952',
                                  '0.60363',
                                  '0.369616',
                                  '0.28265'
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
                                      ]
                  );
@center_rawresults_vector = (
                              'input',
                              1,
                              1,
                              1,
                              1,
                              1,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              3,
                              '74.9062',
                              'FOCE_I',
                              undef,
                              '510.31',
                              '57.67',
                              '10797.653468906055',
                              '1.39797',
                              '4.96413',
                              '4.01979',
                              '-0.262909',
                              '-0.124656',
                              '0.183952',
                              '0.60363',
                              '0.369616',
                              '0.28265',
                              '5e-13',
                              '0.176668',
                              '0.461182',
                              '0.667734',
                              '0.131609',
                              '0.0022974',
                              '0.0753429',
                              '0.18408',
                              '0.0778764',
                              '0.014659',
                              undef,
                              undef,
                              undef,
                              undef,
                              '0.0420175',
                              '0.265058',
                              '0.422124',
                              '0.5161',
                              '0.906451',
                              '1.01602',
                              '1.18166',
                              '1.5032',
                              '3.14737',
                              undef
                            );
@negative_dofv = (
                   0
                 );
@intermediate_raw_results_files = (
                                    'raw_results_final_backward_4.csv'
                                  );
@samples = (
             2000
           );
@resamples = (
               1000
             );
@attempted_samples = (
                       2000
                     );
@successful_samples = (
                        1999
                      );
@actual_resamples = (
                      1000
                    );
@seed_array = (
                2473309,
                0
              );
$with_replacement = 0;
$mceta = 0;
$problems_per_file = 100;
@minimum_ofv = (
                 '10798.696470069339'
               );
$reference_ofv = '10797.653468906055';
$done = 1;
$iteration = 1;
$recenter = 1;
$copy_data = 1;
$boxcox = 1;
$nm_version = 'default';
$model_filename = 'final_backward_4.mod';
$cap_resampling = 1;
$cap_correlation = '0.8';
$adjust_blocks = 0;
$check_cholesky_reparameterization = 0;
$subjects = '75';

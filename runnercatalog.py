0
101 print(f277_a5_list)
102 '''
103 dates = f277_a5_list
104 datetime_objects = [datetime.fromisoformat(date[:-1]) for date in dates]
105
106 # Step 2: Calculate the mean datetime
107 total_timedelta = sum((dt - datetime_objects[0] for dt in datetime_objects), timedelta())
108 mean_timedelta = total_timedelta / (len(datetime_objects) )
109 mean_datetime = datetime_objects[0] + mean_timedelta
110
111 # Step 3: Format the mean datetime back to the desired format
112 mean_datetime_formatted = mean_datetime.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
113
114 print(mean_datetime_formatted)
115 catalog_object = catalog(f277_a5_list[0])
116 catalog_object.concatenate_catalogs(catalog(f277_a5_list[1]), outname='f277w_a5_starcat.fits')
117 catalog_object = catalog('f277w_a5_starcat.fits')
118 catalog_object.concatenate_catalogs(catalog(f277_a5_list[2]), outname='f277w_a5_starcat.fits')
119 catalog_object.concatenate_catalogs(catalog(f277_a5_list[3]), outname='f277w_a5_starcat.fits')
120
121 catalog_object = catalog(f277_b5_list[0])
122 catalog_object.concatenate_catalogs(catalog(f277_b5_list[1]), outname='f277w_b5_starcat.fits')
123 catalog_object = catalog('f277w_b5_starcat.fits')
124 catalog_object.concatenate_catalogs(catalog(f277_b5_list[2]), outname='f277w_b5_starcat.fits')
125 catalog_object.concatenate_catalogs(catalog(f277_b5_list[3]), outname='f277w_b5_starcat.fits')
126
127 catalog_object = catalog(f444_a5_list[0])
128 catalog_object.concatenate_catalogs(catalog(f444_a5_list[1]), outname='f444w_a5_starcat.fits')
129 catalog_object = catalog('f444w_a5_starcat.fits')
130 catalog_object.concatenate_catalogs(catalog(f444_a5_list[2]), outname='f444w_a5_starcat.fits')
131 catalog_object.concatenate_catalogs(catalog(f444_a5_list[3]), outname='f444w_a5_starcat.fits')
132
133 catalog_object = catalog(f444_b5_list[0])
134 catalog_object.concatenate_catalogs(catalog(f444_b5_list[1]), outname='f444w_b5_starcat.fits')
135 catalog_object = catalog('f444w_b5_starcat.fits')
136 catalog_object.concatenate_catalogs(catalog(f444_b5_list[2]), outname='f444w_b5_starcat.fits')
137 catalog_object.concatenate_catalogs(catalog(f444_b5_list[3]), outname='f444w_b5_starcat.fits')
138
139 285     Run PSFEx, creating an output directory if one doesn't already exist.
140 286     Default is to create one directory per exposure.
141 288
142 289     base_name = os.path.basename(image_file)
143 290     outcat_name  = starcat_file.replace('cat.fits', 'psfex_cat.fits')
144 291     psfex_outdir = os.path.join(config['outdir'], \
145 292                     'psfex-output', base_name.split('.')[0])
146 293     if not os.path.isdir(psfex_outdir):
147 294         os.system(f'mkdir -p {psfex_outdir}')
148 295
149 296     psfex_config_arg = '-c ' + os.path.join(config['configdir'],'psfex.config')
150 297     outcat_arg = f'-OUTCAT_NAME {outcat_name}'
151 298     outdir_arg = f'-PSF_DIR {psfex_outdir}'
152 299
153 300     cmd = ' '.join(['psfex',
154 301                 starcat_file, psfex_config_arg, outcat_arg, outdir_arg]
155 302             )
156 303
157 304     print(f'psfex cmd is {cmd}')
158 305     os.system(cmd)
159 306
160
161 directory_path = "/home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs"
162 file_list = glob.glob(directory_path + "/*augmented_starcat.fits")
163
164 print(file_list)
165
166
167 for file in file_list:
168     print(' '.join(['psfex ', file, '-c ' + '/home/eddieberman/research/mcclearygroup/cweb_psf/astro_config/psfex.config', f'-OUTCAT_NAME ' + file.replace('starcat.fits', 'psfex_cat.fits'), f'-PSF_DIR /home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs/']))
169     cmd = ' '.join(['psfex ', file, '-c ' + '/home/eddieberman/research/mcclearygroup/cweb_psf/astro_config/psfex.config', f'-OUTCAT_NAME ' + file.replace('starcat.fits', 'psfex_cat.fits'), f'-PSF_DIR /home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs'])
170     os.system(cmd)
171
172
173  imhead = fits.getheader(image_file, ext=0)
174  24     filter_name = imhead['FILTER']
175  25     date = imhead['DATE']
176  26     detector = imhead['DETECTOR'].replace('LONG', '5')
177  27     roll_angle = fits.getval(image_file, 'ROLL_REF', ext=1)
178  28
179  29     # Define output file name
180  30     if outdir is None:
181  31         outdir = os.path.dirname(image_file)
182  32     image_name = os.path.basename(image_file)
183  33     output_file = os.path.join(outdir,
184  34                     image_name.replace('.fits', '_WebbPSF.fits'))
185  35     # Set oversample scale
186  36     if (filter_name in ['F277W', 'F444W']) and (oversample_lw == True):
187  37         oversample = 2
188  38     else:
189  39         oversample = 1
190  40
191  41     # create WebbPSF instance for given image parameters
192  42     nc = webbpsf.NIRCam()
193  43     nc.options['parity'] = 'odd'
194  44     nc.filter =  filter_name
195  45     nc.oversample = oversample
196  46     nc.detector = detector
197  47
198  48     # Load OSS by date
199  49     nc.load_wss_opd_by_date(date)
200  50
201  51     # Calculate PSF & save to file. Make it big so we don't truncate anything.
202  52     # If crop_psf=False, the distorted image will have different dimensions.
203  53     psf = nc.calc_psf(crop_psf=True, fov_pixels=301)
204  54     print(f"nc.oversample is {nc.oversample}")
205  55     psf.writeto(output_file, overwrite=True)
206  56
207  57     # NOW! We need to rotate the PSF image by the appropriate roll angle
208  58     for i in range(len(psf)):
209  59         image_data = psf[i].data
210  60         rotated_image = rotate(image_data, (360.0-roll_angle))
211  61         psf[i].data = rotated_image
212  62
213  63     psf.writeto(output_file.replace('WebbPSF', 'WebbPSF_rot'), overwrite=True)
214
215
216 def extract_3_numbers(filename):
217     pattern = r'\d{3}'
218     matches = re.findall(pattern, filename)
219     return matches
220 nc = webbpsf.NIRCam()
221
222 for file in file_list:
223     letter_number_codes = re.findall(r'[A-Za-z]\d', file)
224     #print('F'+extract_3_numbers(file)[0]+'W')
225     #print(letter_number_codes[1].replace('a', 'A').replace('b', 'B'))
226     nc.filter = 'F'+extract_3_numbers(file)[0]+'W'
227     nc.detector = 'NRC'+letter_number_codes[1].replace('a', 'A').replace('b', 'B')
228     nc.options['parity'] = 'odd'
229
230     filter_number = extract_3_numbers(file)[0]
231     print(filter_number[0], letter_number_codes[1])
232
233     # Assign the correct list variable based on the filter and detector
234     dates = []
235     if filter_number == '115':
236         if letter_number_codes[1] == 'a1':
237             dates = f115_a1_list
238         elif letter_number_codes[1] == 'a2':
239             dates = f115_a2_list
240         elif letter_number_codes[1] == 'a3':
241             dates = f115_a3_list
242         elif letter_number_codes[1] == 'a4':
243             dates = f115_a4_list
244         elif letter_number_codes[1] == 'b1':
245             dates = f115_b1_list
246         elif letter_number_codes[1] == 'b2':
247             dates = f115_b2_list
248         elif letter_number_codes[1] == 'b3':
249             dates = f115_b3_list
250         elif letter_number_codes[1] == 'b4':
251             dates = f115_b4_list
252     elif filter_number == '150':
253         if letter_number_codes[1] == 'a1':
254             dates = f150_a1_list
255         elif letter_number_codes[1] == 'a2':
256             dates = f150_a2_list
257         elif letter_number_codes[1] == 'a3':
258             dates = f150_a3_list
259         elif letter_number_codes[1] == 'a4':
260             dates = f150_a4_list
261         elif letter_number_codes[1] == 'b1':
262             dates = f150_b1_list
263         elif letter_number_codes[1] == 'b2':
264             dates = f150_b2_list
265         elif letter_number_codes[1] == 'b3':
266             dates = f150_b3_list
267         elif letter_number_codes[1] == 'b4':
268             dates = f150_b4_list
269     elif filter_number == '277':
270         if letter_number_codes[1] == 'a5':
271             dates = f277_a5_list
272         elif letter_number_codes[1] == 'b5':
273             dates = f277_b5_list
274     elif filter_number == '444':
275         if letter_number_codes[1] == 'a5':
276             dates = f444_a5_list
277         elif letter_number_codes[1] == 'b5':
278             dates = f444_b5_list
279
280     print(file)
281     datetime_objects = [datetime.fromisoformat(date[:-1]) for date in dates]
282     total_timedelta = sum((dt - datetime_objects[0] for dt in datetime_objects), timedelta())
283     mean_timedelta = total_timedelta / len(datetime_objects)
284     mean_datetime = datetime_objects[0] + mean_timedelta
285     mean_datetime_formatted = mean_datetime.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
286     print(mean_datetime_formatted)
287     nc.load_wss_opd_by_date(mean_datetime_formatted)
288     nc.oversample = 2
289     psf = nc.calc_psf(crop_psf=True, fov_pixels=301)
290     psf.writeto(file.replace('starcat.fits', 'WebbPSF.fits'), overwrite=True)
291
292
293 for file in file_list:
294     catalog_object = catalog(file)
295     #epsfex_object = epsfex(file.replace('starcat.fits', 'starcat.psf'))
296     print(file.replace('augmented_starcat.fits', 'WebbPSF.fits'))
297     webb_name = file.replace('augmented_starcat.fits', 'WebbPSF.fits')
298     webbpsf_object = webbpsf(webb_name)
299     catalog_object.augment(webbpsf_object, outname=file.replace('augmented_starcat.fits', 'augmented_psf_starcat.fits'))
300 '''
301
302
303
~
~

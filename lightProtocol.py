# File LightProtocol.py
# -*- coding: utf-8 -*-
"""

Copyright (C) 2014-2015  Anna Matuszyńska, Oliver Ebenhöh

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""

class LightProtocol:



    def __init__(self, options={'PFD':0,'protocol':'constant'}):

        PFD = options.setdefault('PFD',0)

        lightFnKey = 'lightintensity'

        proto = options.setdefault('protocol','const')


        if proto.startswith('user'):

            # user defined light function
            # keys in options:
            # 'LightFn': <user_defined_function>
            _userDefinedLightFn = options.setdefault('lightFn',lambda t:0)
            setattr(self,lightFnKey,_userDefinedLightFn)


        elif proto.startswith('const'):

            # constant light
            setattr(self,lightFnKey,lambda t:PFD)


        elif proto == 'ldl':

            # simple light-dark-light protocol
            # required keys in options:
            # 'Toff': <time_at_light_off>, 'Ton': <time_at_light_on>
            options.setdefault('Toff',600);
            options.setdefault('Ton',1200);

            setattr(self,lightFnKey,
                    lambda t: PFD * ((t < options['Toff']) or (t >= options['Ton']))
                    )


        elif proto.startswith('PAM'):

            self.Toff = options.setdefault('Toff',600);
            self.Ton = options.setdefault('Ton',1200);
            self.Tfl = options.setdefault('Tflash',0.8);
            self.dt = options.setdefault('dtFlash',60);
            self.PFDsat = options.setdefault('PFDflash',5000);

            if proto.startswith('PAM_ldl'):
                '''
                LDL with PAM
                '''
                setattr(self,lightFnKey,
                        lambda t: self.PFDsat * ((t % self.dt) < self.Tfl) +
                                  PFD * (((t < self.Toff) or (t > self.Ton)) and (t % self.dt) >= self.Tfl)
                        )

            elif proto.startswith('PAM_dld'):
                '''
                    DLD with PAM
                '''

                setattr(self,lightFnKey,
                        lambda t: self.PFDsat * ((t % self.dt) < self.Tfl) +
                                  PFD * (((t > self.Ton) and (t < self.Toff)) and (t % self.dt) >= self.Tfl)
                        )

        else:
            raise NameError('Protocol \''+options['protocol']+'\' not yet implemented.')

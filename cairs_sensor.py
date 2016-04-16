__author__ = 'matt'
from cairs_coor import Coor


class Sensor(object):
    """
    function that computes log p(S| [R_1 ,... R_n], I)
    must take a signal S and an array of *not transformed* rain
    """

    def __init__(self, log_p, delta_coor=None, domain_extent=None):
        """@:param log_p - methods for computing
           @:param delta_coor - list of Coor
           @:param domain_extent - cube of integration, relative to sensor position
        """
        self.log_p = log_p
        if not delta_coor and not domain_extent:
            coor = Coor(0, 0)
            self.delta_coor = [coor]
            self.domain_extent = coor
        else:
            if delta_coor:
                self.delta_coor = [delta_coor]
                self.domain_extent = Coor(0, 0)
            else:
                self.delta_coor = []
                self.domain_extent = domain_extent


    def __str__(self):
        str1 = 'Sensor(%s,' % self.log_p.__name__
        if self.delta_coor:
            str1 += ' delta_coor:'
            if type(self.delta_coor) is list:
                str1 += '['
                for d in self.delta_coor:
                    str1 += '%s' % str(d)
                    if len(self.delta_coor) > 1:
                        str1 += ', '
                str1 += '],'
            else:
                str1 += '%s' % str(self.delta_coor)
        if self.domain_extent:
            str1 += ' domain_extent:'
            if type(self.domain_extent) is list:
                str1 += '['
                for d in self.domain_extent:
                    str1 += '%s' % str(d)
                    if len(self.domain_extent) > 1:
                        str1 += ', '
                str1 += ']'
            else:
                str1 += '%s' % str(self.domain_extent)

        str1 += ' )'
        '''
        if type(self.delta_coor) is list:
            str += "- log likelihood: %s\n" % self.log_p
            for coor in self.delta_coor:
                str += "- coor %s\n" % coor
                # str+=coor+'\n'

        #if not self.domain_extent.IsZero():
        str += "- integration domain: %s\n" % self.domain_extent
        #else:
        str += "- no integration\n"
        '''
        return str1

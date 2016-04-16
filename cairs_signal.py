from cairs_coor import *


class Signal():
    def __init__(self, signalData, sensor, coor, angle=0):
        """
        :param signal:  signal obtained by sensor
        :param sensor: sensor
        :param coor: coordinates of the sensor
        :param angle: direction in rad, rotations centre is 'coor'
        """
        self.signalData = signalData
        self.sensor = sensor
        self.position = coor
        self.angle = angle

    def __str__(self):

        str = 'Signal!(%s,%s,%s,%s)' % (self.signalData, self.sensor, self.position, self.angle)
        """
        str = "  Signal:\n"
        str += "- position: %s\n" % self.position
        str += "- angle: %s\n" % self.angle
        str += "- measured value: %s \n" % self.signalData
        str += "- signal type: %s\n" % type(self)
        str += "- sensor: %s\n" % self.sensor
        """
        return str

    @staticmethod
    def find_near_signals(loc_pred, signals, delta):

        # Find the extrem time-coordiantes of 'loc_pred'
        st = StaticContext()

        tmin = (datetime.datetime.max - st.REF_TIME).total_seconds()
        tmax = (datetime.datetime.min - st.REF_TIME).total_seconds()
        for loc in loc_pred:  # todo check
            if type(loc) == type(Domain):
                t1 = loc.position.time
                t2 = loc.position.time + loc.extend.time

                tmin = min(tmin, t1, t2)
                tmax = max(tmax, t1, t2)
            else:
                tmin = min(tmin, loc.time)
                tmax = max(tmax, loc.time)

        tmin -= delta
        tmax += delta

        # Find signals that are within (tmin, tmax)
        signals_near = []
        for sig in signals.signals:
            if type(sig.sensor.delta_coor) is list:
                t = [sig.position.time,
                     (sig.position.time + sig.sensor.domain_extent.time),
                     [(sig.position.time + i.time) for i in sig.sensor.delta_coor]]
            else:
                t = [sig.position.time,
                     (sig.position.time + sig.sensor.domain_extent.time),
                ]

            flagYes = False
            for i in t:  # TODO check it
                if type(i) is list:
                    for ii in i:
                        if tmax > ii > tmin:
                            flagYes = True
                else:
                    if tmax > i > tmin:
                        flagYes = True

            if flagYes:
                signals_near.append(sig)
                flagYes = False

        return signals_near


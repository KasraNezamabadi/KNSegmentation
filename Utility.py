import os
import pandas as pd
from scipy import signal
from mat4py import loadmat
import numpy as np
from numpy.linalg import norm
import statistics as stat

path_to_ecg_files = 'Data/ECG'
path_to_ecg_meta = 'Data/ECGMeta.xlsx'
path_to_pla_ann = 'Data/PLAAnnotation'


class WaveBoundaryException(Exception):
    def __init__(self, message):
        super().__init__(message)


class NoOffsetException(WaveBoundaryException):
    def __init__(self, onset_index: int):
        self.onset_index = onset_index
        self.message = 'No Offset for Onset'
        super(NoOffsetException, self).__init__(self.message)


class Util:

    @staticmethod
    def get_offset(onset: int, offset_list: [int], frequency: int):
        # 500 ms is the threshold
        threshold = round(0.5 * frequency)
        closest_index = -1
        for index in range(len(offset_list)):
            if abs(offset_list[index] - onset) < threshold:
                threshold = abs(offset_list[index] - onset)
                closest_index = index
        if closest_index == -1:
            raise NoOffsetException
        return offset_list[closest_index]

    @staticmethod
    def get_min_distance(target: int, input_list: [int]):
        temp = [abs(x - target) for x in input_list]
        return min(temp)

    @staticmethod
    def get_closest_ann(ann_source: int, input_list: [(int, int)]):
        min_diff = 5000
        closest_match = ()
        for item in input_list:
            if min_diff > abs(ann_source - item[0]):
                min_diff = abs(ann_source - item[0])
                closest_match = item

        return closest_match

    @staticmethod
    def get_ecg_list() -> [int]:
        ann_names = [f for f in os.listdir(path_to_pla_ann) if not f.startswith('.')]
        ecg_ids = []
        for ann_name in ann_names:
            ecg_ids.append(int(ann_name.split('_')[0]))
        ecg_ids = list(set(ecg_ids))
        return ecg_ids

    @staticmethod
    def get_lead_name(index: int):
        names = ['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']
        return names[index]


class Loader:

    @staticmethod
    def get_ecg(ecg_id: int, denoise=True) -> [pd.DataFrame, int]:
        meta = pd.read_excel(path_to_ecg_meta)
        try:
            frequency = int(meta.loc[meta['ECG ID'] == ecg_id]['Sample Base'].values[0])
            path = os.path.join(path_to_ecg_files, str(ecg_id) + '.csv')
            ecg = pd.read_csv(filepath_or_buffer=path,
                              header=None, skiprows=1,
                              names=['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'])

            if denoise:
                for lead in ecg:
                    lead_ecg = ecg[lead]
                    filtered = SignalProcessing.filter(sig=lead_ecg, frequency=frequency)
                    ecg[lead] = filtered
            return [ecg, frequency]
        except KeyError:
            assert False, 'Could not find meta data for ECG ' + str(ecg_id)
        except IndexError:
            assert False, 'ECG ' + str(ecg_id) + ' does not have sample frequency.'

    @staticmethod
    def _correct_annotations(ann_list):
        result = []
        for lead_ann in ann_list:
            if type(lead_ann) is not list:
                lead_ann = [lead_ann]
            lead_result = []
            for ann in lead_ann:
                if type(ann) is list:
                    if len(ann) == 0:
                        continue
                    lead_result.append([ann[0] - 1, ann[1] - 1])
                else:
                    lead_result.append(ann-1)
            result.append(lead_result)
        return result

    @staticmethod
    def get_annotations(ecg_id: int) -> [pd.DataFrame]:
        p_ann_name = str(ecg_id) + '_P.mat'
        qrs_ann_name = str(ecg_id) + '_QRS.mat'
        t_ann_name = str(ecg_id) + '_T.mat'
        p_ann = loadmat(os.path.join(path_to_pla_ann, p_ann_name))
        qrs_ann = loadmat(os.path.join(path_to_pla_ann, qrs_ann_name))
        t_ann = loadmat(os.path.join(path_to_pla_ann, t_ann_name))

        p_start = Loader._correct_annotations(p_ann['P_anns']['onset'])
        p_peak = Loader._correct_annotations(p_ann['P_anns']['peak'])
        p_end = Loader._correct_annotations(p_ann['P_anns']['offset'])

        qrs_start = Loader._correct_annotations(qrs_ann['QRS_anns']['onset'])
        qrs_peak = Loader._correct_annotations(qrs_ann['QRS_anns']['peak'])
        qrs_fragmented = Loader._correct_annotations(qrs_ann['QRS_anns']['fragmented'])
        qrs_end = Loader._correct_annotations(qrs_ann['QRS_anns']['offset'])

        t_start = Loader._correct_annotations(t_ann['T_anns']['onset'])
        t_peak = Loader._correct_annotations(t_ann['T_anns']['peak'])
        t_end = Loader._correct_annotations(t_ann['T_anns']['offset'])

        final_ann = []
        for lead in range(12):
            row = [p_start[lead], p_peak[lead], p_end[lead],
                   qrs_start[lead], qrs_peak[lead], qrs_end[lead], qrs_fragmented[lead],
                   t_start[lead], t_peak[lead], t_end[lead]]
            final_ann.append(row)

        final_df = pd.DataFrame(data=final_ann)
        final_df.columns = ['P Start', 'P Peak', 'P End',
                            'QRS Start', 'QRS Peak', 'QRS End', 'QRS Fragmented',
                            'T Start', 'T Peak', 'T End']

        return final_df


class SignalProcessing:
    @staticmethod
    def filter(sig: [int], frequency: int, ):
        nyq = frequency / 2  # Nyquist Frequency
        Wp = [1/nyq, 100/nyq]  # Pass band (Normalised)
        Ws = [0.5 /nyq, 110/nyq]  # Stop band (Normalised)
        order, Wn = signal.cheb2ord(Wp, Ws, gpass=10, gstop=30)
        b, a = signal.cheby2(order, 30, Wn, btype='bandpass', output='ba')
        filtered = signal.filtfilt(b=b, a=a, x=sig)
        return filtered

    def get_significant_points(self, segment: [int], threshold: float):
        significant_points = []
        p_start = (0, segment[0])
        p_end = (len(segment) - 1, segment[-1])
        p_start = np.asarray(p_start)
        p_end = np.asarray(p_end)

        max_distance = 0
        max_index = -1
        for i in range(len(segment)):
            point = (i, segment[i])
            point = np.asarray(point)
            distance = norm(np.cross(p_end - p_start, p_start - point)) / norm(p_end - p_start)
            if distance > max_distance:
                max_distance = distance
                max_index = i

        if max_distance > threshold:
            significant_points.append(max_index)

            left_segment = segment[:max_index]
            right_segment = segment[max_index+1:]

            significant_points_left = self.get_significant_points(segment=left_segment, threshold=threshold)
            significant_points_right = self.get_significant_points(segment=right_segment, threshold=threshold)
            if significant_points_right is not None:
                significant_points_right = [x + max_index for x in significant_points_right]

            if significant_points_left is not None:
                significant_points.extend(significant_points_left)
            if significant_points_right is not None:
                significant_points.extend(significant_points_right)
            return significant_points




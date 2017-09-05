#include <stdio.h>

extern "C" {

#define __STDC_CONSTANT_MACROS

#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswresample/swresample.h>
#include <libavutil/avutil.h>

}

#include <memory>
#include <vector>

namespace {

struct AVFormatCloserAndDeleter {
	void operator() (AVFormatContext *ctx) {
		avformat_close_input(&ctx);
		avformat_free_context(ctx);
	}
};

struct AVCodecContextDeleter {
	void operator() (AVCodecContext *ctx) {
		avcodec_close(ctx);
		av_freep(&ctx);
	}
};

struct SwrContextDeleter {
	void operator() (SwrContext *swr) {
		swr_free(&swr);
	}
};

struct AVPacketDeleter {
	void operator() (AVPacket *pkt) {
		av_free_packet(pkt);
	}
};

#if (LIBAVCODEC_VERSION_MAJOR >= 55)
struct AVFrameDeleter {
	void operator() (AVFrame *frame) {
		av_frame_free(&frame);
	}
};
#endif

struct AVSampleDeleter {
	void operator() (uint8_t *data) {
		av_freep(&data);
	}
};

void convert_samples(SwrContext *swr, int sample_rate, const uint8_t **data, int nb_samples, std::vector<float> *samples)
{
	int max_out_samples = nb_samples + swr_get_delay(swr, sample_rate);
	if (max_out_samples == 0) {
		return;
	}
	uint8_t *output;
	av_samples_alloc(&output, nullptr, 1, max_out_samples, AV_SAMPLE_FMT_FLT, 0);
	std::unique_ptr<uint8_t, AVSampleDeleter> output_deleter(output);

	int out_samples = swr_convert(swr, &output, max_out_samples, data, nb_samples);
	if (out_samples > 0) {
		const float* start = reinterpret_cast<const float *>(output);
		const float* end = start + out_samples;
		samples->insert(samples->end(), start, end);
	}
}

int decode_packet(const char *filename, AVCodecContext *codec_ctx, SwrContext *swr, AVFrame *audio_frame, AVPacket *packet, int *got_frame, std::vector<float> *samples)
{
	*got_frame = 0;
	int len1 = avcodec_decode_audio4(codec_ctx, audio_frame, got_frame, packet);
	if (len1 < 0 || !*got_frame) {
		return len1;
	}

	if (audio_frame->channel_layout != codec_ctx->channel_layout ||
	    audio_frame->sample_rate != codec_ctx->sample_rate) {
		fprintf(stderr, "%s: Channel layout or sample rate changed mid-file\n", filename);
		*got_frame = false;
		return len1;
	}
	convert_samples(swr, codec_ctx->sample_rate, (const uint8_t **)audio_frame->data, audio_frame->nb_samples, samples);
	return len1;
}

}  // namespace

bool read_audio_file(const char *filename, std::vector<float> *samples, int *sample_rate)
{
	av_register_all();

	AVFormatContext *format_ctx = nullptr;
	if (avformat_open_input(&format_ctx, filename, nullptr, nullptr) != 0) {
		fprintf(stderr, "Couldn't open %s\n", filename);
		return false;
	}
	std::unique_ptr<AVFormatContext, AVFormatCloserAndDeleter> format_ctx_closer(format_ctx);

	if (avformat_find_stream_info(format_ctx, nullptr) < 0) {
		fprintf(stderr, "%s: Couldn't find stream information\n", filename);
		return false;
	}

	// Find the first audio stream.
	int audio_stream_index = -1;
	for (unsigned i = 0; i < format_ctx->nb_streams; ++i) {
		if (format_ctx->streams[i]->codec->codec_type == AVMEDIA_TYPE_AUDIO) {
			audio_stream_index = i;
			break;
		}
	}
	if (audio_stream_index == -1) {
		fprintf(stderr, "%s: Couldn't find an audio stream\n", filename);
		return false;
	}

	AVCodec *codec = avcodec_find_decoder(format_ctx->streams[audio_stream_index]->codec->codec_id);
	if (codec == nullptr) {
		fprintf(stderr, "%s: Unsupported codec\n", filename);
		return false;
	}

	AVCodecContext *codec_ctx = avcodec_alloc_context3(codec);
	std::unique_ptr<AVCodecContext, AVCodecContextDeleter> codec_ctx_deleter(codec_ctx);
	if (avcodec_copy_context(codec_ctx, format_ctx->streams[audio_stream_index]->codec) != 0) {
		fprintf(stderr, "%s: Couldn't copy codec context\n", filename);
		return false;
	}

	if (avcodec_open2(codec_ctx, codec, nullptr) < 0) {
		fprintf(stderr, "%s: Couldn't open codec\n", filename);
		return false;
	}

	// Init resampler (to downmix to mono and convert to s16).
	if (codec_ctx->channel_layout == 0) {
		codec_ctx->channel_layout = av_get_default_channel_layout(codec_ctx->channels);
	}
	SwrContext *swr = swr_alloc_set_opts(
		nullptr,
		AV_CH_LAYOUT_MONO, AV_SAMPLE_FMT_FLT, codec_ctx->sample_rate,
		codec_ctx->channel_layout, codec_ctx->sample_fmt, codec_ctx->sample_rate,
		0, nullptr);
	std::unique_ptr<SwrContext, SwrContextDeleter> swr_deleter(swr);
	if (swr_init(swr) < 0) {
		fprintf(stderr, "%s: Couldn't initialize resampler\n", filename);
		return false;
	}

	AVPacket packet;
#if (LIBAVCODEC_VERSION_MAJOR >= 55)
	AVFrame *audio_frame = av_frame_alloc();
	std::unique_ptr<AVFrame, AVFrameDeleter> audio_frame_deleter(audio_frame);
#else
	AVFrame frame_holder {};
	AVFrame *audio_frame = &frame_holder;
#endif
	while (av_read_frame(format_ctx, &packet) >= 0) {
		std::unique_ptr<AVPacket, AVPacketDeleter> av_packet_deleter(&packet);

		if (packet.stream_index != audio_stream_index) {
			continue;
		}

		while (packet.size > 0) {
			int got_frame = 0;
			int len1 = decode_packet(filename, codec_ctx, swr, audio_frame, &packet, &got_frame, samples);
			if (len1 < 0) {
				fprintf(stderr, "%s: Couldn't decode audio\n", filename);
				return false;
			}
			if (!got_frame) {
				break;
			}
			packet.data += len1;
			packet.size -= len1;
		}
	}

	// Flush any delayed data from the end.
	packet.data = nullptr;
	packet.size = 0;
	int got_frame = 0;
	do {
		int len1 = decode_packet(filename, codec_ctx, swr, audio_frame, &packet, &got_frame, samples);
		if (len1 < 0) {
			fprintf(stderr, "%s: Couldn't decode audio\n", filename);
			return false;
		}
	} while (got_frame);

	// Convert any leftover samples from the converter.
	convert_samples(swr, codec_ctx->sample_rate, nullptr, 0, samples);

	*sample_rate = codec_ctx->sample_rate;

	return true;
}
